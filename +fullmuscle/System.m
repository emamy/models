classdef System < models.muscle.System;
% System: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-16
%
% @new{0,7,dw,2014-09-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=private)
        num_motoneuron_dof;
        num_sarco_dof;
        num_spindle_dof;
        num_all_dof;
        num_extra_dof;
        off_moto;
        off_sarco;
        off_spindle;
        
        input_motoneuron_link_idx;
        moto_input_noise_factors;
        sarco_output_idx;
        
        % The offset for the sarcomere signals at t=0. Used to create
        % alpha(X,0) = 0
        %
        % See also: assembleX0
        sarco_mech_signal_offset
        
        Motoneuron;
        Spindle;
        Sarcomere;
        HasSpindle;
    end
    
    properties(SetAccess=private)
        nfibres;
    end
    
    methods
        function this = System(model)
            this = this@models.muscle.System(model);
            this.f = models.fullmuscle.Dynamics(this);
            
            mc = model.Config;
            opts = mc.Options;
            
            % First row is neumann input
            this.Inputs{1,1} = mc.getAlphaRamp(30,1);
            % Second row is external mean current input
%             this.Inputs{2,1} = @(t)1;
            
%             % Setup noise input
%             ng = models.motoneuron.NoiseGenerator;
%             ng.DisableNoise = ~options.Noise;
%             this.noiseGen = ng;
%             this.Inputs{1} = @ng.getInput;
            
            if opts.SarcoVersion == 1
                this.Sarcomere = models.motorunit.SarcomereOriginal;
            else
                this.Sarcomere = models.motorunit.SarcomereNeumann;
            end
            this.num_sarco_dof = this.Sarcomere.Dims;
            
            this.Motoneuron = models.motoneuron.Motoneuron;
            this.num_motoneuron_dof = this.Motoneuron.Dims;
            
            % Compile information for plotting
            this.Plotter = models.fullmuscle.MusclePlotter(this);
            
            % Constraint nonlinearity
            this.g = models.fullmuscle.ExtendedConstraint(this);
        end
        
        function configUpdated(this)
            mc = this.Model.Config;
            ft = mc.FibreTypes;
            nf = length(ft);
            this.nfibres = nf;
            
            %% Configure motoneurons
            this.Motoneuron.setType(ft);
            
            %% Configure sarcomeres
            this.Sarcomere.setType(ft);
            
            %% Configure Spindle (if set)
            hassp = ~isempty(mc.SpindlePositions);
            this.HasSpindle = hassp;
            this.Spindle = [];
            if hassp
                this.Spindle = models.fullmuscle.Spindle;
            end
            
            configUpdated@models.muscle.System(this);
        end
        
        function setConfig(this, mu, inputidx)
            setConfig@models.muscle.System(this, mu, inputidx);
            
            if ~isempty(inputidx)
                % Create an input substitute that uses the true external
                % function and creates the effective noisy signal from it
                maxcurrents = this.Motoneuron.getMaxMeanCurrents;
                
                % Get noise from motoneuron class
                no = this.Motoneuron.TypeNoise;
                bno = this.Motoneuron.BaseNoise;
                
                % First row is neumann input
                uneum = this.Inputs{1,inputidx};
                
                % Second row is external mean current input
                uext = this.Inputs{2,inputidx};
                
                ustr = '@(t)[mu(3)*uneum(t); bno(round(t)+1); ';
                for k=1:this.nfibres
                    rowfun = sprintf('no(%d,round(t)+1)*min(%g,mu(4)*uext(t)); ', k, maxcurrents(k));
                    ustr = [ustr rowfun];%#ok
                end
                ustr = [ustr ']'];
                this.u = eval(ustr);
            end
        end
        
        function uvwall = includeDirichletValues(this, t, uvw)
            uvwall_mech = includeDirichletValues@models.muscle.System(...
                this, t, uvw(1:this.off_moto,:));
            uvwall = [uvwall_mech; uvw(this.off_moto+1:end,:)];
        end
        
    end
    
    methods(Access=protected)
        
        function updateDimensions(this, mc)
            updateDimensions@models.muscle.System(this, mc);
            
            dm = this.Motoneuron.Dims;
            dsa = this.Sarcomere.Dims;
            
            % Motoneurons are beginning after mechanics
            this.off_moto = this.NumTotalDofs;
            this.num_motoneuron_dof = dm*this.nfibres;

            % Sarcomeres are beginning after motoneurons
            this.off_sarco = this.off_moto + this.num_motoneuron_dof;
            this.num_sarco_dof = dsa*this.nfibres;
            
            % Spindles are beginning after sarcomeres
            this.off_spindle = this.off_sarco + this.num_sarco_dof;
            this.num_spindle_dof = 0;
            if this.HasSpindle
                ds = this.Spindle.Dims;
                this.num_spindle_dof = ds*this.nfibres;
            end
            this.num_all_dof = this.off_spindle + this.num_spindle_dof;
            
            % Get the positions where the input signal is mapped to the
            % motoneurons
            this.input_motoneuron_link_idx = this.off_moto + (2:dm:dm*this.nfibres);
            
            this.sarco_output_idx = this.off_sarco + (53:dsa:dsa*this.nfibres);
            
            % Set all those extra dofs as algebraic dofs
            this.num_extra_dof = this.num_motoneuron_dof...
                + this.num_sarco_dof ...
                + this.num_spindle_dof;
            this.NumAlgebraicDofs = this.NumAlgebraicDofs ...
                + this.num_extra_dof;
            % Re-compute the total dof number
            this.NumTotalDofs = this.NumStateDofs ...
                + this.NumDerivativeDofs + this.NumAlgebraicDofs;
        end
        
        function ad_ic_as_x0 = getAlgebraicDofsInitialConditions(this)
            ad_ic_as_x0 = getAlgebraicDofsInitialConditions@models.muscle.System(this);
            
            dm = this.Motoneuron.Dims;
            dsa = this.Sarcomere.Dims;
            mc = this.Model.Config;
            opt = mc.Options;
            
            % Actual constraint dofs (=g operator from mechanics)
            off = this.NumAlgebraicDofs - this.num_extra_dof;
            
            % Load dynamic/affine x0 coefficients for moto/sarco system
            % from file
            meta = metaclass(this);
            s = load(fullfile(fileparts(which(meta.Name)),...
                    sprintf('x0coeff%d.mat',opt.SarcoVersion)));
            x0_motorunit = dscomponents.AffineInitialValue;
            m = size(s.coeff,1);
            for k=1:m
                x0_motorunit.addMatrix(sprintf('polyval([%s],mu(1))',...
                    sprintf('%.14g ',s.coeff(k,:))),full(sparse(k,1,1,m,1)));
            end
            % Evaluate the dynamic ICs for the selected fibre types and use
            % them as x0
            ft = mc.FibreTypes;
            smoff = zeros(this.nfibres,1);
            for k=1:this.nfibres
                x0ms = x0_motorunit.evaluate(ft(k));
                % add moto
                ad_ic_as_x0(off + dm*(k-1) + (1:dm)) = x0ms(1:dm);
                % add sarco
                off = off + this.num_motoneuron_dof;
                ad_ic_as_x0(off + dsa*(k-1) + (1:dsa)) = x0ms((dm+1):end);
                smoff(k) = x0ms(dm+53);
                if this.HasSpindle
                    ds = this.Spindle.Dims;
                    off = off + this.num_sarco_dof;
                    % add spindle
                    ad_ic_as_x0(off + ds*(k-1) + (1:ds)) = this.Spindle.y0;
                end
            end
            this.sarco_mech_signal_offset = smoff;
        end
        
%         function x0 = assembleX0(this)
%             x0 = zeros(this.NumTotalDofs,1);
%             % Get muscle x0 - fills state space until motoneuron offset
%             x0(1:this.NumStateDofs) = assembleX0@models.muscle.System(this);
%         end
        
%         function Baff = assembleB(this)
%             Baff = dscomponents.AffLinInputConv;
%             
%             %% Add neumann forces B into larger affine matrix in first
%             % column
%             % Assemble B matrix for neumann conditions from models.muscle.System
%             Baff_neumann = assembleB@models.muscle.System(this);
%             if ~isempty(Baff_neumann)
%                 B = sparse(this.NumTotalDofs,this.nfibres+2);
%                 B(1:this.off_moto,1) = Baff_neumann.getMatrix(1);
%                 Baff.addMatrix(Baff_neumann.funStr{1},B);
%             end
%             
% %             %% Add motoneuron external input signal
% %             i = this.input_motoneuron_link_idx;
% %             s = this.Motoneuron.FibreTypeNoiseFactors;
% %             % Initialize with base noise entries in second column
% %             B = sparse(i,ones(size(i))+1,s,this.num_all_dof,this.nfibres+2);
% %             % Add single contribution for each fibre type thereafter
% %             for k=1:this.nfibres
% %                 B(i(k),k+2) = s(k);%#ok
% %             end
% %             Baff.addMatrix('1',B);
%         end
    end
    
end