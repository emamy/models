classdef System < models.BaseFirstOrderSystem
% MotoSystem: The global dynamical system used within the MotoModel
%
%
% @author Daniel Wirtz @date 2012-11-22
%
% @new{0,7,dw,2013-07-05} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(Transient)
        % debug output
        d = [];
    end
    
    properties (Constant)        
        dm = 6; % Dimension of motoneuron part
        
        % Membrane capacitance. Different values for different fibre types, 
        % due to different action potential propagation speeds
        % C_m is computed parameter-dependent.
        % These constants are for both slow and fast muscle models and are also used in the
        % first entry of the sarcomere constants computed in
        % models.muscle.FibreDynamics.initSarcoConst @type double
        C_m_slow = 0.58;
        C_m_fast = 1;
    end
    
    properties(SetAccess=private)
        % The NoiseGenerator to obtain the inputs u(t)      
        noiseGen;
    end
    
    methods
        function this = System(model)
            % Call superclass constructor
            this = this@models.BaseFirstOrderSystem(model);
            
            this.addParam('fibre_type', 0, 'Range', [0 1],'Desired',30);
            this.addParam('mean_current_factor', 3, 'Range', [0 9],'Desired',30);
            % Extra parameter only for NoiseScaling experiment!
%             this.addParam('noise_scaling', 1, 'Range', [.8 2],...
%                 'Spacing','log','Desired',30);
            
            this.NumStateDofs = 6;
            
            % Setup noise input
            ng = models.motoneuron.NoiseGenerator;
            this.noiseGen = ng;
            this.Inputs{1} = @ng.getInput;
            
            %% Set system components
            % Core nonlinearity
            this.f = models.motoneuron.Dynamics(this);
            
            % Linear input B for motoneuron
            % input conversion matrix, depends on fibre type. Has only one
            % entry in second row.
            B = dscomponents.AffLinInputConv;
            
            % Mean input mapping - scales with µ_2 as mean current factor
            B.addMatrix('(mu(2,:))./(pi*(exp(log(100)*mu(1,:))*3.55e-05 + 77.5e-4).^2)',...
                sparse(2,1,1,6,2));
            % Noise input mapping - scales with fractional root of µ_2 
            % over [0,9] as factor
            % See models.motoneuron.experiments.NoiseScaling!
            B.addMatrix('(9*((mu(2,:)/9).^1.24635))./(pi*(exp(log(100)*mu(1,:))*3.55e-05 + 77.5e-4).^2)',...
                sparse(2,2,1,6,2));
            this.B = B;
            
            % Dynamic initial values - same as computed for motorunit
            mc = ?models.motorunit.Shorten;
            s = load(fullfile(fileparts(which(mc.Name)),'x0coeff1.mat'));
            x0 = dscomponents.AffineInitialValue;
            for k=1:6
                x0.addMatrix(sprintf('polyval([%s],mu(1))',...
                    sprintf('%g ',s.coeff(k,:))),full(sparse(k,1,1,6,1)));
            end
            this.x0 = x0;
            %this.x0 = dscomponents.ConstInitialValue(zeros(6,1));
            
            this.updateDimensions;
        end
        
        function prepareSimulation(this, mu, inputidx)
            % Limits the mean input current to a maximum value depending on
            % the fibre type, so that the frequency of 60Hz is not
            % exceeded.
            %
            % See also: models.motoneuron.ParamDomainDetection upperlimit_poly
            
            % Limit mean current depending on fibre type
            if this.Model.FibreTypeDepMaxMeanCurrent
                mu(2) = this.f.moto.checkMeanCurrent(mu(2),mu(1));
            end
            prepareSimulation@models.BaseFirstOrderSystem(this, mu, inputidx);
        end
        
        function setConfig(this, mu, inputidx)
            setConfig@models.BaseFirstOrderSystem(this, mu, inputidx);
            this.noiseGen.setFibreType(mu(1));
            this.MaxTimestep = this.Model.dt*1000;            
        end
        
        function pm = plot(this, t, y, pm)
            if nargin < 4
                n = 2;
                if ~isempty(this.mu)
                    n = 3;
                end
                pm = PlotManager(false,1,n);
                pm.LeaveOpen = true;
            end
            
            % Input
            h = pm.nextPlot('input','Input currents (base + indep + mean)','time','value');
            plot(h,t,this.noiseGen.getInput(t));
            
            % Effective input
            if ~isempty(this.mu)
                B = this.B.compose([],this.mu);
                noise = B(2,:)*this.noiseGen.getInput(t);
                h = pm.nextPlot('eff_input',...
                    sprintf('Effective input current, noise mean=%g',mean(noise)),'time','value');
                plot(h,t,noise);
            end
            
            % Firing rate
            [peaks, cov_isi,cov_dr,freq] = this.Model.analyze(t, y);
            h = pm.nextPlot('moto',...
                sprintf('Motoneuron V_s\nmu=[%g %g]\nCoV(ISI)=%g, CoV(DR)=%g, Hz=%g',...
                this.mu,cov_isi,cov_dr,freq),'time','value');
            plot(h,t,y(2,:),'b',t(peaks),y(2,peaks),'rx');
            
            if nargin < 4
                pm.done;
            end
        end
    end    
end
