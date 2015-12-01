classdef Model < models.BaseFullModel
    % Model: Electromyography model
    %
    % Model representing the second equation in the bidomain equations.
    % `\nabla\cdot[(\sigma_i+\sigma_e)\nabla\phi_e] + \nabla\cdot[\sigma_i\nabla\phi_e] = 0`
    %
    % @author Timm Strecker @date 2014-03-24
    %
    % @change{0,8,dw,2015-09-21} Imported into +models package.
    %
    % @new{0,7,ts,2014-03-24} Added this class.
    %
    
    properties(Constant)
        DataDir = fullfile(fileparts(mfilename('fullpath')),'data');
    end
    
    properties
        % Seed that can be used by random number generator instances in order to enable result
        % reproduction.
        %
        % @type int @default 1
        RandSeed = 1;
    end
    
    properties(Dependent)
        % Position (index) of the neuromuscular junction of each fibre. By
        % default randomly anywhere on the fibres.
        % matrix<integer> of dimension dim(2) \times zdim_m
        neuronpos;
    end
    
    properties(SetAccess=private)
        Options;
    end
    
    methods
        function this = Model(varargin)
            % Creates a new Electromyography model
            %
            % Geo: vector with 4 entries representing length,
            % width, muscle height and height of the fat/skin layer (in cm). If
            % there are only 3 entries, skinheight is set to zero
            % Dimension of the discretization in all 3 dimensions.
            
            rs = RandStream('mt19937ar','Seed',12345);
            
            i = inputParser;
            i.addParameter('Dim',[40;20;13],@(v)numel(v)==3);
            i.addParameter('Geo',0.1*[40;20;10;3],@(v)numel(v)==3 || numel(v)==4);
            i.addParameter('MUTypes',[0 rs.rand(1,3) 1]); % 5 types by default
            i.addParameter('Shapes','precomp',...
                @(v)any(strcmp(v,{'actual','precomp','rosen','full'})));
            i.addParameter('FiringTimes','precomp',@(v)any(strcmp(v,{'actual','precomp'})));
            % The sarcomere implementation version for computation of
            % shorten's action potential shapes
            i.addParameter('SarcoVersion',1);
            i.addParameter('DynamicAmplitudes',true,@(v)islogical(v) && isscalar(v));
            i.addParameter('DynamicPropagationSpeed',true,@(v)islogical(v) && isscalar(v));
            i.parse(varargin{:});
            opts = i.Results;
            
            this.Options = opts;
            % Make sure we have column vectors
            geo = reshape(opts.Geo,[],1);
            opts.Dim = reshape(opts.Dim,[],1);
            if numel(geo) == 3
                geo = [geo; 0];
            end
            opts.Geo = geo;
            
            this.Name = sprintf('EMG model geometry = [%g,%g,%g,%g] discretization = [%g,%g,%g]',geo, opts.Dim);
            this.System = models.emg.System(this, opts);
            this.SaveTag = sprintf('EMG_geom[%g,%g,%g,%g]_discr[%g,%g,%g]', geo, opts.Dim);
            this.Data = data.ModelData(this);
            
            this.isStatic = true;
            this.DefaultMu = 5;
            
            d = approx.DEIM(this.System);
            d.MaxOrder = 2000;
            this.Approx = d;
            s = spacereduction.PODReducer;
            s.Value=1e-5;
            s.Mode='eps';
            this.SpaceReducer=s;
            
            this.EnableTrajectoryCaching = true;
            
            this.T = 100;
            this.dt = 1;
        end
        
        function distributeMotorUnits(this, varargin)
            % Distributes number many motor units over muscle-crosssection.
            % All fibres of each MU lie within a circle. By default, the
            % circles are distributed approximately uniformly over the
            % muscle-crosssection such that every fibre can be assigned to
            % at least one MU. The optional arguments centers and radii can
            % be used to define custom circles.
            %
            % Parameters:
            % number: The number of motor units
            % @type integer
            % centers: centers of the circles @type 2 \times number matrix<double>
            % radii: radii of the circles @type vector<double> of length
            % number or scalar
            this.System.f.distributeMotorUnits(varargin{:});
        end
        
        function initNeuronposGaussian(this, varargin)
            % Initializes neuron position by sampling from a Gaussian
            % normal distribution standard deviation sigma and mean mu. 
            % By default, the mean value is in the middle of the muscle 
            % fibres and the standard deviation equals 0.3cm.
            this.System.f.initNeuronposGaussian(varargin{:});
        end
        
        %% Model reduction methods
        function [reduced,time] = buildReducedModel(this, target_dim, DEIM_order)
            tic;
            if nargin < 3
                DEIM_order = find(this.Approx.SingularValues > 1e-7, 1, 'last');
                if nargin < 2
                    %target_dim = size(this.Data.V,2);
                    target_dim = find(this.SpaceReducer.SingularValues > 1e-4, 1, 'last');
                end
            end
            
            reduced = buildReducedModel@models.BaseFullModel(this, target_dim);
            reduced.System.f.Order = DEIM_order;
            time = toc;
        end
        
        function W = computeWfromV(this, V)%#ok
            % Idea 1: (Timm) use W = A*V
%             W = this.System.A.evaluate(1,0)*V;
            
            % Idea 2: Simply use zeros for the last row - it's a side
            % condition which is satisfied for the reduced basis anyways
            % (by definition as the values are generated from trajectories
            % that satisfy zero mean at each point, hence does the POD)
            W = [V; zeros(1,size(V,2))];
            
            o = general.Orthonormalizer;
            W = o.orthonormalize(W);
        end
        
        function off4_genApproximationTrainData(this)
            % Use W=V for off4_genApproximationTrainData, otherwise there 
            % is an error in data.ApproxTrainData.computeFrom
            % (because V and W don't have the same number of rows).
            oldW = this.Data.W;
            this.Data.W = this.Data.V;
            off4_genApproximationTrainData@models.BaseFullModel(this);
            this.Data.W = oldW;
        end
        
        function off5_computeApproximation(this)
            A = KerMor.App;
            host = A.getHost;
            
            % On LEAD: Write ApproxTrainData.fxi in a matrix - this makes
            % POD of fxi much faster
            if length(host) >= 4 && strcmp(host(1:4), 'lead')
                start = tic;
                atd = this.Data.ApproxTrainData;
                this.Approx.computeDEIM(this.System.f, atd.fxi.toMemoryMatrix);
                this.OfflinePhaseTimes(5) = toc(start);
            else
                off5_computeApproximation@models.BaseFullModel(this);
            end
        end
    end
    
    %% plotting methods 
    % illustrating the motor unit configuration
    methods  
        function pm = plot(this, t, y, pm)
            if nargin < 4
                pm = PlotManager;
                pm.LeaveOpen = true;
            end
            sys = this.System;
            xaxis = 0:sys.h(1):sys.Geo(1);
            yaxis = 0:sys.h(2):sys.Geo(2);
            ax = pm.nextPlot('emg','','width [cm]','length [cm]');
            zlabel(ax,'\phi_i [mV]');
            pm.done;
            view(53,28);
            %axis([0 this.Geo(2) 0 this.Geo(1) min(hlp(:)) max(hlp(:))]);
            maxval = ceil(max(y(:)));
            minval = floor(min(y(:)));
            axis([0 sys.Geo(2) 0 sys.Geo(1) minval maxval]);
            caxis([minval maxval]);
            % daspect([1 1 1]);
            daspect([1 1 (maxval-minval)]);
            hold(ax,'on');
            for idx = 1:length(t)
                cla(ax);
                phi = reshape(y(:,idx),sys.dim(1),sys.dim(2),[]);
                % surf(ax, yaxis, xaxis, 10*phi(:,:,end),'FaceColor','interp','EdgeColor','interp');  % pcolor?
                surf(ax, yaxis, xaxis, phi(:,:,end),'FaceColor','interp','EdgeColor','interp');  % pcolor?
                title(ax,sprintf('Surface EMG at time %gms',t(idx)));
                colorbar;
                colormap jet;
                pause(.005);
                if ~ishandle(ax)
                    break;
                end
            end
            if nargin < 4
                pm.done;
            end
        end
        
        function plotMuscleConfiguration(this)
            this.plotMUcrosssection;
            mg = this.System.Geo;
            fprintf('The muscle is %gcm long, %gcm wide and %gcm high.\n',mg(1),mg(2),mg(3));
            if (mg(4) ~= 0)
                fprintf('There is a %gcm high skin/ fat layer on top of the msucle.\n',mg(4));
            end
            h = this.System.h;
            fprintf('The discretization width is h = [%g, %g, %g]cm in x-, y- and z-direction.\n',h(1),h(2),h(3));
            fprintf('The muscle fibres are grouped into %g motor units.\n',length(this.System.f.MUTypes));
            FT = this.Options.FiringTimes;
            if strcmpi(FT,'precomp')
                disp('The motoneurons'' firing times were precomputed by a motoneuron model.');
            elseif strcmp(FT,'actual')
                disp('The motoneurons'' firing times are determined by simulation of a motoneuron model.');
            else
                disp('Custom motoneuron firing times are used.');
            end
        end
        
        function plotMU3d(this, numbers)
            % plots all fibres of one motor unit over the muscle.
            % number: index of the motor unit
            %rhs = this.System.f;
            sys = this.System;
            numstr = '';
            figure
            hold on
            if nargin < 2
                numbers = 1:min(length(this.System.f.MUTypes), 7);
            end
            ls = LineSpecIterator(length(numbers));
            for num = numbers
                pos=[];
                for i=1:sys.dim(2)
                    for j=1:sys.zdim_m
                        if sys.f.MUGeoTypeIdx(i,j) == num
                            pos = [pos,[i;j]];%#ok
                        end
                    end
                end
                for i = 1:size(pos,2)
                    X = [0,sys.Geo(1)];
                    Y = [(pos(1,i)-1)*sys.h(2),(pos(1,i)-1)*sys.h(2)];
                    Z = [(pos(2,i)-1)*sys.h(3), (pos(2,i)-1)*sys.h(3)];
                    plot3(X, Y, Z, 'Color',ls.nextColor,'LineWidth',2)
                    plot3((this.neuronpos(pos(1,i),pos(2,i))-1)*sys.h(1),Y(1),Z(1),'xk','LineWidth',2)
                end
                numstr = [numstr, sprintf(', #%g',num)];%#ok
            end
            grid on
            
            title(['Muscle Fibres and Neuromuscular Junctions of Motor Unit(s) ',numstr]);
            xlabel('x-direction [cm]');
            ylabel('y-direction [cm]');
            zlabel('z-direction [cm]');
            axis([0 sys.Geo(1) 0 sys.Geo(2) ...
                0 sys.Geo(3)+sys.Geo(4)])
            view(-44,12);
            daspect([1,1,1]);
            hold off
        end
        
        function plotMUcrosssection(this, numbers)
            % plots the motor unit distribution over muscle cross section.
            % Musclefibres of the same motor unit have the same color and
            % marker.
            
            rhs = this.System.f;
            nmu = length(rhs.MUTypes);
            if nargin < 2  % plot all motor units
                numbers = 1:nmu;
            end
            
            sys = this.System;
            geo = sys.Geo;
            pos=cell(1,nmu);
            for i=1:sys.dim(2)
                for j=1:sys.zdim_m
                    pos{rhs.MUGeoTypeIdx(i,j)}=[pos{rhs.MUGeoTypeIdx(i,j)},[i;j]];
                end
            end
            figure
            % frames for body region
            if geo ~= 0
                X = [0, 0, geo(2), geo(2), 0];
                Y = [geo(3), geo(3)+geo(4), ...
                    geo(3)+ geo(4), geo(3), geo(3)];
                patch(X', Y', [0;0;0;0;0], 'FaceColor',[0.9 0.9 .9])
            end
            hold on
            
            ls = LineSpecIterator(length(numbers),1000);
            for num = numbers
                if isempty(pos{num})
                    continue;
                end
                col = ls.nextColor;
                plot((pos{num}(1,:)-1)*sys.h(2),(pos{num}(2,:)-1)*sys.h(3),...
                    '.', 'Marker', ls.nextMarkerStyle, 'MarkerEdgeColor', col, 'LineWidth', 1);
                if ~isempty(rhs.MUCenters) && nmu > 1
                    plotCircle(rhs.MUCenters(:,num),rhs.MURadii(num));
                end
            end
            title('Motor Unit Distribution over Muscle Cross Section');
            xlabel('y-direction [cm]');
            ylabel('z-direction [cm]');
            axis([0,geo(2), 0,geo(3)+geo(4)])
            daspect([1 1 1]);
            hold off
            
            function plotCircle(cent, rad)
                dphi = 0.01;
                phi = 0:dphi:2*pi;
                x = cent(1) + rad*cos(phi);
                y = cent(2) + rad*sin(phi);
                if geo(4) > 0
                    idx = logical(y <= geo(3)*1.01);
                    x = x(idx);
                    y = y(idx);
                end
                plot(x,y,'-','Color',col);
            end
        end
        
        function h = plot1D(this, t, y, point)
            % plots the EMG signal at point over time
            geo = this.System.Geo;
            if nargin < 4
                point = geo(1:3)/2;
            end
            if numel(point) ~= 3
                error('Input argument point must have 3 entries.\n');
            end
            if any([point(1)>geo(1), point(2)>geo(2), point(3)>geo(3)+geo(4)]) || any(point < 0)
                error('Input argument point must lie inside muscle.\n');
            end
            
            idx(1) = round(point(1)/this.System.h(1)) + 1; % ./ fails if point is row and h is column
            idx(2) = round(point(2)/this.System.h(2)) + 1;
            idx(3) = round(point(3)/this.System.h(3)) + 1;
            dim = this.System.dim;
            idx = idx(1) + dim(1)*(idx(2)-1) + dim(1)*dim(2)*(idx(3)-1);
            
            h = figure;
            plot(t, y(idx,:), t, abs(y(idx,:)),'LineWidth',2);
            title(sprintf('EMG signal at point (%g,%g,%g)',point(1),point(2),point(3)));
            xlabel('time [ms]');
            yl = '\phi_e';
            if point(3) > geo(3)
                yl = '\phi_o';
            end
            ylabel([yl, ' [mV]']);
            legend([yl, '(t)'], ['|',yl,'(t)|']);
        end
        
        function plotAPShapes(this)
            rhs = this.System.f;
            nmu = length(rhs.MUTypes);
            pm = PlotManager(false,2,2);
            pm.LeaveOpen = true;
            for k=1:nmu                
                ft = rhs.MUTypes(k);
                ax = pm.nextPlot(sprintf('apshape%d_type%g',k,ft),...
                    sprintf('ActionPotential shape %d: fibre_type %g',k,ft),...
                    't [ms]','Vm [mV]');
                sh = rhs.APShapes{k};
                t = (0:(length(sh)-1))*this.dt;
                plot(ax,t,sh);
            end
            pm.done;
        end
        
        function pm = plotPrecompAPShapes(this)
            sv = this.Options.SarcoVersion;
            % Else: "precomp" is set, so load & interpolate shapes
            % for current time.
            datafile = fullfile(fileparts(mfilename('fullpath')),'data',...
                sprintf('ShapeData_v%d.mat',sv));
            s = load(datafile);
            pm = PlotManager;
            pm.LeaveOpen = true;
            ax = pm.nextPlot('precomp_apshapes',...
                sprintf('Precomputed AP shapes for %d parameters',length(s.Times)),...
                't [ms]','potential [mV]');
            hold(ax,'on');
            for k=1:length(s.Times)
                plot(ax,s.Times{k},s.Shapes{k});
            end
        end
        
        function pm = plotVm(this, mean_current)
            if nargin < 2
                mean_current = this.DefaultMu;
            end
            sys = this.System;
            xdim = sys.dim(1);
            sys.prepareSimulation(mean_current,[]);
            types = this.Options.MUTypes;
            nmu = length(types);
            [T,X] = meshgrid(this.Times,(-xdim+1:xdim-1)*sys.h(1));
            pm = PlotManager(false,1,2);
            pm.LeaveOpen = true;
            for k=1:nmu
                Vm = sys.f.computeMUActivation(this.Times,k);
                % Use only the "right" side (symmetric)
                %Vm = Vm(xdim:end,:);
                ax = pm.nextPlot('Vm',...
                    sprintf('Vm signal over space and time for mu=%g and fibre type %g',...
                    mean_current,types(k)),...
                    't [ms]','fibre location [cm]');
                surf(T,X,Vm,'Parent',ax,'EdgeColor','interp');
            end
            pm.done;
        end
    end
    
    methods  % getters and setters: interface to System and RHS
        function value = get.neuronpos(this)
            value = this.System.f.neuronpos;
        end
        
        function set.neuronpos(this, value)
            if ~all(size(value) == [this.System.dim(2) this.System.zdim_m])
                error('Wrong dimension.')
            end
            this.System.f.neuronpos = value;
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'models.emg.Model')
                sobj = this;
                this = models.emg.Model;
                this.RandSeed = sobj.RandSeed;
                this = loadobj@models.BaseFullModel(this, sobj);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end
end