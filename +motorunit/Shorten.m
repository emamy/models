classdef Shorten < models.BaseFullModel
    % Shorten: Model for a muscle motor unit composed of motoneuron
    % and a sarcomere.
    %
    % The global time unit for this model is milliseconds [ms].
    % This model has been copied from the MuscleFibreModel in order to
    % obtain a current snapshot without changing the multi-sarcomere
    % adoptions done for the whole fibre model.
    %
    % This model is intended for research regarding the
    % kernel-approximation of the motor unit w.r.t. the activation and
    % fibre type relation to the force development.
    %
    % @author Daniel Wirtz @date 2014-01-16
    %
    % @new{0,7,dw,2014-01-16} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Dependent)
        UseNoise;
    end
    
    methods
        function this = Shorten(varargin)
            % Creates a new motor unit model
            
            i = inputParser;
            i.addParameter('SarcoVersion',1);
            i.addParameter('DynamicIC',true);
            i.addParameter('SPM',false);
            i.addParameter('OutputScaling',true);
            i.addParameter('Noise',true);
            i.parse(varargin{:});
            options = i.Results;
            
            this.dt = .1;
            if options.SPM
                this.T = 2000; % [ms]
            else
                this.T = 150; % [ms]
            end
            
            this.dt = 0.1; % [ms]
            
            this.SaveTag = sprintf('motorunit_shorten_sp%d_dynic%d',...
                options.SPM,options.DynamicIC);
            this.Data = data.ModelData(this);
            this.Data.useFileTrajectoryData;
            
            this.Name = sprintf('Motor unit model: Single peak mode: %d, Dynamic Initial Conditions: %d',...
                options.SPM,options.DynamicIC);
            this.System = models.motorunit.System(this, options);
            this.TrainingInputs = 1;
            this.EnableTrajectoryCaching = true;
            
            this.ODESolver = solvers.MLWrapper(@ode15s);
            
            % Set parameter domain already
            s = sampling.RandomSampler;
            s.Domain = models.motoneuron.ParamDomain;
            this.Sampler = s;
            
            this.DefaultMu = [.1; 3];
            this.DefaultInput = 1;
        end
        
        function plotOutputForceScaling(this, x)
            if nargin < 2
                x = 0:.01:1;
            end
            plot(x,polyval(this.System.ForceOutputScalingPolyCoeff,x));
            title('Force scaling curve for different fibre types');
            xlabel('Fibre type parameter');
            ylabel('Peak force for single excitation');
        end
        
        function pm = plotState(this, t, x, pm)
            if nargin < 4
                pm = PlotManager(false,3,1);
                pm.LeaveOpen = true;
            end
            h = pm.nextPlot('moto','Motoneuron V_s','time','value');
            plot(h,t,x(2,:));
            h = pm.nextPlot('sarco','Linked sarcomere: V_s','time','V_s');
            plot(h,t,x(this.System.dm+1,:));
            h = pm.nextPlot('sarco',sprintf('Linked sarcomere: A_2\nMu=[%s]',num2str(this.System.mu')),'time','A_2');
            plot(h,t,x(this.System.dm+53,:));
            
            if nargin < 4
                pm.done;
            end
        end
        
        function pm = plot(this, t, y, pm)
            if nargin < 4
                pm = PlotManager(false,2,1);
                pm.LeaveOpen = true;
            end
            h = pm.nextPlot('sarco','Linked sarcomere: V_s','time','V_s');
            plot(h,t,y(1,:));
            h = pm.nextPlot('sarco',sprintf('Linked sarcomere: A_2\nMu=[%s]',num2str(this.System.mu')),'time','A_2');
            plot(h,t,y(2,:));
            
            if nargin < 4
                pm.done;
            end
        end
        
        function [apshape, times, ct] = getActionPotentialShape(this, fibre_type, basemV)
            % Computes the action potential shape for the current Shorten
            % model. Always detects the first peak.
            %
            % We recommend using DynamicIC=true and SPM=true for speed.
            %
            % Parameters:
            % fibre_type: The desired fibre type within [0,1] @type double
            % @default this.DefaultMu(1)
            % basemV: The base potential from which to measure the shape.
            % The shape's start will be the first time-step where the
            % signal is greater than basemV and the end will be the first
            % time-step the signal will be below basemV.
            % @type double @default -80[mV]
            %
            % Return values:
            % apshape: The shape of the action potential with values as-is
            % from the Shorten model. @type rowvec<double>
            % times: The corresponding times over which the shape is
            % computed. @type rowvec<double>
            % ct: The computation time in seconds @type double
            if nargin < 3
                basemV = -80;
                if nargin < 2
                    fibre_type = this.DefaultMu(1);
                end
            end
            if this.T < 40
                warning('T < 40ms! Action potential shape might be uncorrect!');
            end
            % Always run with max activation to immediately have a peak
            [t,y,ct] = this.simulate([fibre_type;9],1);
            thlp = tic;
            % Criteria via absolute value is bad for slow-twitch fibres.
            startidx = find((y(1,:) > basemV),1,'first');
            diffstartidx = find((diff(y(1,:)) > .01),1,'first');
            maxlowms = .1;
            if (diffstartidx - startidx)*this.dt > maxlowms
                startidx = max(1,diffstartidx-round(maxlowms/this.dt));
            end
            len = find((y(1,startidx:end) < basemV),1,'first');
            pos = startidx + (1:len);
            apshape = y(1,pos);
            times = t(pos)-min(t(pos));
            ct = ct + toc(thlp);
        end
    end
    
    methods
        function value = get.UseNoise(this)
            value = ~this.System.noiseGen.DisableNoise;
        end
        
        function set.UseNoise(this, value)
            this.System.noiseGen.DisableNoise = ~value;
        end
    end
    
    methods(Static)
        function res = test_Shorten
            res = 1;
            m = models.motorunit.Shorten;
            [t,y] = m.simulate;
            m.plot(t,y);
            
            m = models.motorunit.Shorten('SarcoVersion',2);
            [t,y] = m.simulate;
            m.plot(t,y);
            
            m.System.f.LinkSarcoMoto = false;
            [t,y] = m.simulate;
            m.plot(t,y);
            
            m = models.motorunit.Shorten('SarcoVersion',1,...
              'DynamicIC',false,'SPM',true);
            m.T = 150;
            [t,y] = m.simulate;
            m.plot(t,y);
            
            for sv=1:2
                m = models.motorunit.Shorten('SarcoVersion',sv,...
                      'DynamicIC',false,'SPM',true);
                m.T = 150;
                m.simulate;
                
                for mu=[0 1]
                    m = models.motorunit.Shorten('SarcoVersion',sv,'SPM',true,...
                        'OutputScaling',true);
                    m.T = 150;
                    [~,y] = m.simulate([mu; 9]);
                    % Need zero force at beginning
                    if ~isequal(y(2,1),0)
                        fprintf('Initial force not equal to zero!');
                        res = false;
                    end
                    % The scaling needs to be such that a single peak has
                    % force equal to one. Allow 1% difference
                    % See the +experiments/SarcoScaling script.
                    if max(y(2,:)-1)>.01
                        fprintf('More than 1%% error on force scaling!');
                        res = false;
                    end
                end
            end
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'models.motorunit.Shorten')
                sobj = this;
                this = models.motorunit.Shorten;
                this = loadobj@models.BaseFullModel(this, sobj);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end 
end