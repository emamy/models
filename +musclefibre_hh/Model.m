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
       
    properties
        Dimension;
        A_m;
        C_m;
        Sigma_eff;
        Dx;
        L;
   end
    
    properties(Access=private)
        fDim;
    end
    
    methods
        function this = Model(dim, am, cm, sigma_eff, dx, t, dt, l)            
%             name = sprintf('Muscle fibre model (%d cells)',options.N);
            this = this@models.BaseFullModel();
            if nargin < 1
                dim = 2048;
            end
            % wird bei Burgers im Hauptfile festgelegt
            this.T = t; % [ms]
            this.dt = dt; % [ms]
            this.A_m = am;
            this.C_m = cm;
            this.Sigma_eff = sigma_eff;
            this.Dx = dx;
            this.L = l;
            this.System = models.musclefibre_hh.System(this);
            this.ODESolver = solvers.SemiImplicitEuler(this);
            this.Dimension = dim;
%             this.SaveTag = sprintf('musclefibre_%dcells',options.N);
            % this.Data = data.ModelData(this);
            this.System.MaxTimestep = this.dt;
            this.System.MaxTimestep = this.dt;
            this.SaveTag = sprintf('musclefibre_1D_d%d',dim);
            this.Name = sprintf('1D-%dd hodkin-huxley equation',dim);
%             if options.N > 1000
%                 this.Data.useFileTrajectoryData;
%             end    
            this.SpaceReducer = spacereduction.PODGreedy;
%             this.SpaceReducer.Eps = 1e-9;
            

            a = approx.DEIM(this.System);
            a.MaxOrder = 80;
            this.Approx = a;
            
            this.ErrorEstimator = error.DEIMEstimator;
            % mache ich sp�ter
%             this.Sampler.Domain = models.motoneuron.ParamDomain;
%             
%             this.DefaultMu = [.1; 3];
%             this.DefaultInput = 1;
        end
        
        function set.Dimension(this, value)
            this.fDim = value;
            this.System.newDim;
            this.Data.SimCache.clearTrajectories;
        end
        
        function dim = get.Dimension(this)
            dim = this.fDim;
        end
        
        

        
        
        function plot(this, t, y, pm_ax)
            if nargin < 4
                pm_ax = PlotManager;
                pm_ax.LeaveOpen = true;
            end
            nt = length(t);
            y  = [zeros(nt,1) y' zeros(nt,1)]; % add boundaries
            xx = linspace(0, size(y,2)-1, size(y,2));
            
            if ~ishandle(pm_ax)
                pm_ax = pm_ax.nextPlot('musclefibre',sprintf('%s: dim=%d',this.Name,this.fDim),'x','t');
            end
            surfc(pm_ax,xx,t,y);
            shading interp;
            zlabel('y');
            rotate3d on;
        end
    end
    

    
%     methods(Static, Access=protected)
%         function this = loadobj(this)
%             if ~isa(this, 'models.musclefibre_hh.Model')
%                 sobj = this;
%                 this = models.musclefibre.Model(sobj.System.N);
%                 this.RandSeed = sobj.RandSeed;
%                 this = loadobj@models.BaseFullModel(this, sobj);
%             else
%                 this = loadobj@models.BaseFullModel(this);
%             end
%         end
%     end
end