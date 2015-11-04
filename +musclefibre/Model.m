classdef Model < models.BaseFullModel
    % Model: Model for a muscle fibre model composed of motoneuron, spindle and
    % sarcomere array.
    %
    % The global time unit for this model is milliseconds [ms].
    %
    %
    % ToDos:
    % - Jacobimatrix-implementierung: evaluatePartialDerivatives. Dann könnten die eingebauten
    %  solver ggf schneller lösen, und der fehlerschätzer könnte implementiert werden
    % - Parameterbereich einschränken: nicht [0,1] sondern vielleicht erstmal [0,.2].
    % - Feststellen, welches der insgesamt zu betrachtende Simulationszeitraum sein sollte (abh.
    % vom aktuellen parameterbereich)
    % - Versuch, größere simulation für N>=400 oder so anzuwerfen
    % - Motoneuron untersuchen und parameterbereich auf sinnvolle größen einschränken! dazu eigenes
    % moto-model bauen und parameterbereich abtasten, dann parameterdomänen-geometrie einbauen.
    %
    % @author Daniel Wirtz @date 2012-11-22
    %
    % @new{0,8,dw,2015-09-15} Imported into +models package.
    %
    % @new{0,7,dw,2012-11-22} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(Constant)
        dxDefault = 10^-2.5;
    end
    
    properties
        % Seed that can be used by random number generator instances in order to enable result
        % reproduction.
        % @type int @default 1
        RandSeed = 1;
    end
    
    properties(SetAccess=private)
        Options;
    end
    
    methods
        function this = Model(varargin)
            % Creates a new muscle fibre model
            %
            % Parameters:
            % N: number of sarcomeres @default 100
            
            i = inputParser;
            i.addParameter('SarcoVersion',1);
            i.addParameter('DynamicIC',true);
            i.addParameter('SPM',false);
            i.addParameter('OutputScaling',true);
            i.addParameter('Spindle',false);
            i.addParameter('N',100,@(n)n>2);
            i.addParameter('Noise',true,@(v)islogical(v));
            i.addParameter('JunctionN',1,@(n)isposintscalar(n));
            % For default dx, see the
            % experiments.PropagationSpeed_SpatialConvergence* scripts
            i.addParameter('dx',models.musclefibre.Model.dxDefault,...
                @(v)isscalar(v));
            i.parse(varargin{:});
            options = i.Results;
            
            if options.JunctionN > options.N
                error('The junction N must be smaller or equal to N');
            end
            
            name = sprintf('Muscle fibre model (%d cells)',options.N);
            this = this@models.BaseFullModel(name);
            this.Options = options;
            this.System = models.musclefibre.System(this, options);
            this.TrainingInputs = 1;
            this.SaveTag = sprintf('musclefibre_%dcells',options.N);
            this.Data = data.ModelData(this);
            if options.N > 1000
                this.Data.useFileTrajectoryData;
            end
            
            this.T = 400; % [ms]
            % DO NOT INCREASE! Peaks from Motoneuron are not resolved
            % visually correct if larger timesteps are used, leading to
            % confusion...
            this.dt = .1; % [ms]
            
            this.ODESolver = solvers.MLWrapper(@ode15s);
            
            this.Sampler.Domain = models.motoneuron.ParamDomain;
            
            this.DefaultMu = [.1; 3];
            this.DefaultInput = 1;
        end
        
        function pm = plot(this, t, y, pm)
            if nargin < 4
                pm = PlotManager(false,1,2);
                pm.LeaveOpen = true;
            end
            h = pm.nextPlot('moto','Motoneuron V_m','time','value');
            plot(h,t,y(1,:));
            o = this.Options;
            mu = this.System.mu;
            if isempty(mu)
                mu = this.DefaultMu;
            end
            tit = sprintf('Sarcomere (Version=%d) A_2, mu=%g\nN=%d, SPM=%d, DynIC=%d, OS=%d',...
                o.SarcoVersion,mu(1),o.N,o.SPM,o.DynamicIC,o.OutputScaling);
            h = pm.nextPlot('sarco',tit,'time [ms]','fibre position [cm]');
            dx = this.System.dx*(1:(size(y,1)-1));
            [DX,T] = meshgrid(t,dx);
            surf(DX,T,y(2:end,:),'EdgeColor','none','FaceColor','interp','Parent',h);
            hold(h,'on');
            Jpos = o.JunctionN+1;
            o = ones(size(t));
            plot3(h,t,dx(Jpos)*o,0*o,'r');
            if nargin < 4
                pm.done;
            end    
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'models.musclefibre.Model')
                sobj = this;
                this = models.musclefibre.Model(sobj.System.N);
                this.RandSeed = sobj.RandSeed;
                this = loadobj@models.BaseFullModel(this, sobj);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end
end