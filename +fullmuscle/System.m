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
        Motoneuron;
        Spindle;
        Sarcomere;
        HasSpindle;
    end
    
    properties(SetAccess=private)
        nfibres;
    end
    
    methods
        function this = System(model, mc)
            this = this@models.muscle.System(model);
            % Here we have a force argument for muscle contraction!
            %this.HasForceArgument = true;
            
            this.f = models.fullmuscle.Dynamics(this);
            
            opts = mc.Options;            
            if opts.SarcoVersion == 1
                this.Sarcomere = models.motorunit.SarcomereOriginal;
            else
                this.Sarcomere = models.motorunit.SarcomereNeumann;
            end
            
            this.Motoneuron = models.motoneuron.Motoneuron;
            
            % Constraint nonlinearity - same as for normal muscle model
            this.g = models.muscle.Constraint(this);
            
            % Add the motorunits and spindle dynamics as standard first
            % order dynamics
            this.FO = models.fullmuscle.FirstOrderDynamics(this);
        end
        
        function configUpdated(this)
            mc = this.Model.Config;
            if ~isempty(mc)
                ft = mc.FibreTypes;
                this.nfibres = length(ft);

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

                this.FO.configUpdated;
            end
        end
        
%         function setConfig(this, mu, inputidx)
%             setConfig@models.muscle.System(this, mu, inputidx);
%             
%             if ~isempty(inputidx)
%                 % Create an input substitute that uses the true external
%                 % function and creates the effective noisy signal from it
%                 maxcurrents = this.Motoneuron.getMaxMeanCurrents;
%                 
%                 % Get all type-dependent noises from motoneuron class
%                 no = this.Motoneuron.Noise;%#ok
%                 
%                 % First row is neumann input
%                 uneum = this.Inputs{1,inputidx};%#ok
%                 
%                 % Second row is external mean current input
%                 uext = this.Inputs{2,inputidx};%#ok
%                 
%                 % Neumann input as first dimension
%                 % Motoneuron base mean as second
%                 % Motoneuron type-dep noises as third
%                 ng = models.motoneuron.NoiseGenerator;
%                 ustr = sprintf('@(t)[mu(3)*uneum(t); %g*ones(size(t)); ',ng.AP);
%                 for k=1:this.nfibres
%                     rowfun = sprintf('no(%d,round(t)+1)*min(%g,mu(4)*uext(t)); ', k, maxcurrents(k));
%                     ustr = [ustr rowfun];%#ok
%                 end
%                 ustr = [ustr ']'];
%                 this.u = eval(ustr);
%             end
%         end
        
    end
    
    methods(Access=protected)
        
        function updateDimensions(this, mc)
            % Recompute the first order dimensions (number of fibre types
            % could have changed)
            fo = this.FO;
            fo.updateDimensions(mc);
            this.NumFirstOrderDofs = fo.fDim;
            
            updateDimensions@models.muscle.System(this, mc);
        end
        
        function x0 = getFOx0(this, ~)
            % We put the x0 assembly into the first order dynamics class as
            % all the relevant quantities are stored there
            x0 = this.FO.x0;
        end
    end
    
end