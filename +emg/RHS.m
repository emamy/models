classdef RHS < dscomponents.ACompEvalCoreFun
    % FibreDynamics: Class for right-hand side B*u(t) of the EMG model
    %
    % @author Timm Strecker @date 2014-04-08
    %
    % @change{0,8,dw,2015-09-21} Imported into +models package.
    %
    % @new{0,7,ts,2014-04-08} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(SetAccess = private)
        B; % the B matrix (discrete generalized Laplace with conductivity sigma_i)
        
        % velocity of action potentials propagating along the muscle
        % fibres in m/s (=> 10cm/ms)
        APvelocity = 2;
        
        % cell array with precomputed action potential shapes for each
        % fibre type
        APShapes;
        % Minimum and range values for each shape
        APShapes_misc;
        
        MUCenters;
        MURadii;
        MUTypes;
        MUFiringTimes;
        
        % distribution of motor units over cross section
        %
        % @default ones(dim(2),zdim_m)
        MUGeoTypeIdx;
    end
    
    properties 
        % Position (index) of the neuromuscular junction of each fibre. By
        % default randomly anywhere on the fibres.
        % matrix<integer> of dimension dim(2) \times zdim_m
        % see also: initNeuronposGaussian
        neuronpos;
    end
    
    properties(SetAccess=private)
        % The randstream instance
        rs;
        
        dims;
        
        % The selected options
        options;
        
        % motoneuron model for computation of firing times at FiringTimesMode 'simulate'
        motomodel;
        last_mean_current;
        
        % Shorten model for action potential shape computation
        shapemodel = [];
        last_dt = [];
        
        % Fields for dynamic amplitude & propagation speed data
        amp_kexp = [];
        ps_kexp = [];
        amp_xiscale = [];
        ps_xiscale = [];
        upperlimit_poly = [];
        
        fullshapes = false;
        dvm;
    end
    
    methods
        function this = RHS(sys)
            this = this@dscomponents.ACompEvalCoreFun(sys);
            this.rs = RandStream('mt19937ar','Seed',24247247);
        end
        
        function init(this, opts)
            this.options = opts;
            this.MUTypes = opts.MUTypes;
            
            this.xDim = prod(opts.Dim);
            this.fDim = prod(opts.Dim)+1;
            sys = this.System;
            this.dims = [sys.dim(1);sys.dim(2);sys.zdim_m];
            
            this.JSparsityPattern = sparse(this.fDim, this.xDim);
            
            this.assembleB;
            
            % Assigns this.MU**
            this.distributeMotorUnits(opts.MUTypes);
            
            this.initDynamicAmpPS(opts);
            this.initAPShapes(opts);
            this.initFiringTimes(opts);
            this.initNeuroJunctions;
            if strcmp(opts.Shapes,'full')
                this.fullshapes = true;
                this.dvm = models.emg.DetailedVm(sys.Geo(1),...
                    sys.dim(1),opts.SarcoVersion);
            end
        end
        
        function copy = clone(this)
            % Create new instance ToDo
            copy = models.emg.RHS(this.System);
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            % Copy local properties
            copy.motomodel = this.motomodel; % dont clone the motoneuron model
            copy.last_mean_current = this.last_mean_current;
            copy.shapemodel = this.shapemodel;
            copy.last_dt = this.last_dt;
            copy.dims = this.dims;
            copy.APShapes = this.APShapes;
            copy.APvelocity = this.APvelocity;
            copy.MUGeoTypeIdx = this.MUGeoTypeIdx;
            copy.neuronpos = this.neuronpos;
            copy.MUCenters = this.MUCenters;
            copy.MURadii = this.MURadii;
            copy.MUTypes = this.MUTypes;
            copy.MUFiringTimes = this.MUFiringTimes;
        end
        
        function proj = project(this, ~, W)
            proj = this.clone;
            proj.B = W'*this.B;
        end
    end
    
    
    %% methods that compute the right-hand side
    methods
        function dy = evaluate(this, ~, t)
            % Evaluates the nonlinear core function at given time. Can evaluate vectorized arguments.
            %
            % Parameters:
            % y: @type matrix<double>, not used
            % t: current time @type rowvec<double>
            % mu: mean input current @type rowvec<double>
            
            dy = this.B * this.getVm(t);
            %if ~isempty(this.W)
            %    dy = this.W'*dy;
            %end
        end
        
        function dy = evaluateMulti(this, ~, t, mu)
            % evaluate can already handle vectorized arguments
            mu = mu(1:this.System.ParamCount,:);
            if all(all(repmat(mu(:,1),1,size(mu,2)-1) == mu(:,2:end)))
                this.prepareSimulation(mu(:,1));
                dy = this.evaluate([], t);
            else
                dy = zeros(this.fDim, length(t));
                for k = 1:length(t)
                    this.prepareSimulation(mu(:,k));
                    dy(:,k) = this.evaluate([], t(k));
                end
            end
        end
        
        function u = getVm(this, t)
            % Computes the "source term" V_m at given time points t.
            % see also: computeMUActivation
            
            d = this.dims;
            u = zeros(prod(d),length(t));
            % Assemble signal for each fibre type
            for muidx = 1:length(this.MUTypes)
                % Compute the activation of that type over length and time
                Vm = this.computeMUActivation(t, muidx);
                
                % Find cross-section positions of current fibre type
                mutype_fibre_idx = find(this.MUGeoTypeIdx == muidx);
                
                % Get positions of neuromuscular junctions
                neuro_junction_pos = this.neuronpos(mutype_fibre_idx);
                
                % Compute the spatial index vectors for each fibre of that
                % type, offset to the junction index
                % The Vm_type is 2*d(1)-1, with the "origin" being index
                % d(1). So, e.g. for junction position 1 we need 40:79 and
                % for position 40 we need 1:40.
                activation_indices = bsxfun(@plus,d(1)-neuro_junction_pos, 1:d(1));
                
                % Compute the positions within u for the fibres of current
                % type incl. correct offset
                u_offsets = bsxfun(@plus,(mutype_fibre_idx-1)*d(1), 1:d(1));
                
                % Magic!
                u(u_offsets(:),:) = Vm(activation_indices(:),:);
            end
        end
        
        function Vm = computeMUActivation(this, t, muidx)
            % Computes the "source term" V_m for all motor units at given
            % time points t. For each MU, V_m is computed over a
            % representative fibre that is twice as long as the actual muscle and
            % has its neuromuscular junction in the middle (=index d(1))
            tlen = length(t);
            xdim = this.dims(1);
            
            ftype = this.MUTypes(muidx);
            mc = this.mu;
            % Restrict mean current to upper limit values!
            % HACK: -1 as the upperlimit poly is okay w.r.t to the
            % resulting frequencies, but the ParamDomain from Timm
            % Strecker is more restrictive. As a consequence, the
            % usable parameters for the learned expansion only make
            % sense within the ParamDomain, which is approximately
            % bounded above by polyval-1.
            mc = min(polyval(this.upperlimit_poly,ftype),mc);
            mu = [ftype; mc];
            if this.fullshapes
                sig = this.dvm.computeSignal(t, mu);
                Vm = zeros(2*xdim-1,tlen);
                Vm(xdim:end,:) = sig;
                Vm(xdim:-1:1,:) = sig;
            else
                o = this.options;    
                
                base = this.APShapes_misc(1,muidx);
                shape = this.APShapes{muidx};
                shapelen = length(shape);

                % Get firing times and remove those that fire beyond
                % the max requested times
                ft = this.MUFiringTimes{muidx}-.4;
                ft(ft > max(t)) = [];

                % Dynamic amplitudes
                if o.DynamicAmplitudes
                    % S is the scaling of mu_1,mu_2,t arguments to the learned
                    % expansion - the time has a different order of magnitude
                    % and is hence scaled so that all arguments are of equal
                    % importance to the expansion (more stable learning)
                    s = this.amp_xiscale;
                    % Arguments: fibre_type, mean_current, time
                    % Important here is to take the speed and amplitudes
                    % for the time instant the junction fire signal was
                    % received, as propagating potentials would have decreasing
                    % amplitudes over time even though the signal would be the
                    % same
                    xi = [repmat(mu./s(1:2),1,length(ft))
                          ft/s(3)];
                    % Get learned amplitudes!
                    amp = this.amp_kexp.evaluate(xi);
                    amp_fac = (amp-base)/this.APShapes_misc(2,muidx);
                else
                    amp_fac = ones(size(ft));
                end
                % Dynamic propagation speeds
                if o.DynamicPropagationSpeed
                    % See o.DynamicAmplitudes case for comments
                    s = this.ps_xiscale;
                    xi = [repmat([ftype; mc]./s(1:2),1,length(ft))
                          ft/s(3)];
                    % Get learned velocities!
                    v = this.ps_kexp.evaluate(xi)/10;
                else
                    v = this.APvelocity/10*ones(size(ft));
                end
                sys = this.System;
                dx = sys.h(1);
                xpos = (0:xdim-1)*dx;
                % Initialize the Vm matrix with double spatial size, so that
                % the signals can be placed in both directions on the fibre
                Vm = base*ones(2*xdim-1,tlen);
                for fidx = 1:length(ft)
                    % Scale shape by current ft amplitude
                    shape_amp = (shape-base)*amp_fac(fidx)+base;
                    % Get the amount of time passed until the signal reaches
                    % each sarcomere (including the local one at 0)
                    shifted_ft = xpos / v(fidx);
                    % Get pairwise difference and find absolute minima =
                    % matching locations
                    [~,shifted_tidx] = min(abs(bsxfun(@minus,t-ft(fidx),shifted_ft')),[],2);
                    for xidx = 1:xdim
                        % Detect how many shape time-steps can be inserted
                        % (might go over the time-index)
                        shapesteps = min(shapelen,tlen-shifted_tidx(xidx));
                        % Set the shapes from center (=xdim) to both sides
                        Vm(xdim+xidx-1,shifted_tidx(xidx)+(1:shapesteps)) = ...
                            shape_amp(1:shapesteps);
                        Vm(xdim-xidx+1,shifted_tidx(xidx)+(shapesteps:-1:1)) = ...
                            shape_amp(shapesteps:-1:1);
                    end
                end
            end
        end
        
        function assembleB(this)
            % Assembles the B matrix of the right-hand side Bu(t).
            % The upper rows of B equal the discrete generalized Laplace \nabla\sigma_i\nabla.
            % These rows correspond to points in the muscle tissue, where there are muscle fibres.
            % All other rows are all-zero. These correspond to the points in
            % the skin layer (no V_m). The last zero-row is due to the
            % zero-mean condition.
            musclegrid = general.geometry.RectGrid3D(this.dims(1),this.dims(2),this.dims(3));
            [iu, ju, su] = find(MatUtils.generalizedLaplacemat3D(this.System.h,musclegrid,this.System.sigmai));
            this.B = sparse(iu, ju, su, this.fDim, prod(this.dims));
        end
        
        function evaluateCoreFun(~, ~, ~, ~)
            error('evaluate is overridden directly.');
        end
        
        function J = getStateJacobian(this, ~, ~, ~)
            % The Jacobian is empty because RHS uses no state space arguments.
            J = sparse(this.fDim,this.xDim);
        end
    end
    
    
    %% methods preparing a simulation
    methods
        function prepareSimulation(this, mu)
            prepareSimulation@dscomponents.ACoreFun(this, mu);
            % For a full fibre simulation we need evaluations of the
            % different musclefibre types
            if ~this.fullshapes
                % We need to compute the current
                % action potential shapes, possibly for each current model
                % time-step dt
                this.updateAPShapes;

                % Update the firing times
                this.updateFiringTimes(mu);
            end
        end
    end
    
    methods(Access=private)
        
        function initAPShapes(this, opts)
            % For actual shape computation use a Shorten model.
            if strcmp(opts.Shapes,'actual')
                m = models.motorunit.Shorten(...
                    'SarcoVersion',opts.SarcoVersion,...
                    'SPM',true,'DynamicIC',true);
                m.T = 50;
                this.shapemodel = m;
            end
        end
        
        function updateAPShapes(this)
            % Choose same step size for precomputed shapes
            fm = this.System.Model;
            cur_dt = fm.dt;
            % Only compute if not yet done or step size changed
            if isempty(this.last_dt) || ~isequal(this.last_dt,cur_dt)
                types = this.MUTypes;
                ntypes = length(types);
                this.APShapes = cell(1,length(types));
                
                % If shapemodel is set, we have the "actual" mode.
                switch this.options.Shapes
                    case 'actual'
                    m = this.shapemodel;
                    pi = ProcessIndicator('Computing %d action potential shapes...',...
                        ntypes,false,ntypes);
                    m.dt = cur_dt;
                    for n = 1:ntypes
                        % Get shape from internal Shorten model using full
                        % activation (for early peak)
                        this.APShapes{n} = m.getActionPotentialShape(types(n));
                        %this.APShapes{n} = interp1(times, shape, ltimes, 'pchip', shape(1));
                        pi.step;
                    end
                case 'precomp'
                    sv = this.options.SarcoVersion;
                    % Else: "precomp" is set, so load & interpolate shapes
                    % for current time.
                    datafile = fullfile(fm.DataDir,...
                        sprintf('ShapeData_v%d.mat',sv));
                    if exist(datafile,'file') ~= 2
                        error('No precomputed data available. Run the %s script!',...
                            fullfile(fm.DataDir,'ActionPotentialShape.m'));
                    end
                    d = load(datafile);
                    pi = ProcessIndicator('Interpolating %d action potential shapes...',...
                        ntypes,false,ntypes);
                    for n = 1:ntypes
                        % Get index of closest matching fibre type (could
                        % be exact match!)
                        [~, idx] = min(abs(d.fibretypes - types(n)));
                        shape = d.Shapes{idx};
                        times = d.Times{idx};
                        % Get absolute time span of shape
                        tspan = times(end)-times(1);
                        % Get that resolved by current time-step
                        ltimes = 0:cur_dt:tspan;
                        % Interpolate on local time grid
                        this.APShapes{n} = interp1(times, shape, ltimes, 'pchip', shape(1));
                        pi.step;
                    end
                    pi.stop;
                case 'rosen'
                    % This is the rosenfalck' ap shape function
                    r = general.functions.Rosenfalck;
                    r.basemV = -80;
                    rosenfalck = r.getFunction;
                    % Average ap shape goes up to 8ms
                    z = 0:cur_dt:8;
                    for n = 1:ntypes
                        % Scale the rosenfalck function by two as the
                        % function is defined over 0-15ms (too long)
                        this.APShapes{n} = rosenfalck(z);
                    end
                end
                
                % Postprocessing (for dynamic shape amplitudes)
                this.APShapes_misc = zeros(2,ntypes);
                for n = 1:ntypes
                    sh = this.APShapes{n};
                    % Compute offset and reference amplitude
                    base = min(sh);
                    this.APShapes_misc(1,n) = base;
                    this.APShapes_misc(2,n) = max(sh-base);
                end
                this.last_dt = cur_dt;
            end
        end
        
        function initFiringTimes(this, opts)
            % init motoneuron model
            if strcmpi(opts.FiringTimes,'actual')
                moto = models.motoneuron.Model(true);
                moto.dt = 0.1;
                moto.EnableTrajectoryCaching = true;
                this.motomodel = moto;
            end
        end
        
        function updateFiringTimes(this, mean_current)
            % Computes the motoneuron firing times.
            %
            % Either pre-computed or for the current fibre-type selection
            % and parameter (=mean input current)
            fm = this.System.Model;
            nmu = length(this.MUTypes);
            if ~isequal(this.last_mean_current,mean_current)
                this.MUFiringTimes = cell(1,nmu);
                m = this.motomodel;
                if ~isempty(m)    
                    m.T = fm.T;
                    pi = ProcessIndicator('Computing firing times for %d motoneurons...',...
                        nmu,false,nmu);
                    for idx = 1:nmu  
                        type = this.MUTypes(idx);
                        [t,y] = m.simulate([type; mean_current], 1);
                        APtimes = (y(2,2:end) > 40).*(y(2,1:end-1) <= 40);
                        this.MUFiringTimes{idx} = t(logical(APtimes));
                        pi.step;
                    end
                    pi.stop;
                else
                    datafile = fullfile(fm.DataDir,'FiringTimes.mat');
                    if exist(datafile,'file') ~= 2
                        error('No precomputed data available. Run the %s script!',...
                            fullfile(fm.DataDir,'FiringTimes.m'));
                    end
                    d = load(datafile);
                    for idx = 1:nmu
                        type = repmat(this.MUTypes(idx),1,size(d.mus,2));
                        [~,pos] = min(Norm.L2(type - d.mus(1,:)));
                        this.MUFiringTimes{idx} = d.Times{pos};
                    end
                end
                this.last_mean_current = mean_current;
            end
        end
        
        function distributeMotorUnits(this, types, centers, radii)
            % Distributes 'num' many motor units over muscle-crosssection.
            % All fibres of each MU lie within a circle. By default, the
            % circles are distributed approximately uniformly over the
            % muscle-crosssection such that every fibre can be assigned to
            % at least one MU. The optional arguments centers and radii can
            % be used to define custom circles.
            %
            % Parameters:
            % num: The number of motor units @type integer @default 1
            % centers: centers of the circles @type rowvec<double>
            % radii: radii of the circles @type rowvec<double>
            
            if nargin < 2
                num = 4;
                % Have slow/fast twitch plus intermediate ones
                types = [0 this.rs.rand(1,num-2) 1];
            else
                num = length(types);
            end
            geo = this.System.Geo;
            if nargin <= 3
                % Randomly compute centers
                centers = bsxfun(@times,geo(2:3),this.rs.rand(2,num));
            else
                if ~all(size(centers) == [2 num])
                    error('centers must be matrix of dimension 2 x %g', num);
                end
            end
            if nargin < 4
                % Daniel: radii = [norm(geo) norm(geo)*(this.rs.rand(1,number-1)*.5 + .1)];
                radii=norm(geo(2:3))*sqrt(exp(types * log(100))/100);
            else
                if isscalar(radii)
                    radii = radii*ones(1,num);
                elseif length(radii) ~= num
                    error('radii must have %g entries or be scalar.', num);
                end
            end
            
            % Store for later use (plots etc)
            this.MUCenters = centers;
            this.MURadii = radii;
            
            % Assign fibres to motor units
            % Choose the maximum random value as fibre for a MU. The chance
            % outside the radius around the circle is set to zero.
            h = this.System.h;
            d = this.dims;
            sel = zeros(d(2),d(3),num);
            [X,Y] = meshgrid(0:h(3):geo(3),0:h(2):geo(2));
            for i=1:num
                rand = this.rs.rand(d(2),d(3));
                in_circle = sqrt((X-centers(1,i)).^2 + (Y-centers(2,i)).^2) < radii(i);
                rand(~in_circle) = 0;
                sel(:,:,i) = rand;
            end
            [~, this.MUGeoTypeIdx] = max(sel,[],3);
        end
        
        function initNeuroJunctions(this, sigma, mu)
            % Initializes neuron position by sampling from a Gaussian
            % normal distribution standard deviation sigma and mean mu.
            % By default, the mean value is in the middle of the muscle
            % fibres and the standard deviation equals 0.3cm.
            sys = this.System;
            if nargin < 3
                mu = sys.Geo(1)/2;
                if nargin < 2
                    sigma = min(0.3, mu/2);
                end
            end
            r = mu + sigma*this.rs.randn(sys.dim(2), sys.zdim_m);
            while any(any(r < 0)) || any(any(r > sys.Geo(1)))
                idx = logical((r < 0) + (r > sys.Geo(1)));
                rnew = mu + sigma*this.rs.randn(sys.dim(2), sys.zdim_m);
                r(idx) = rnew(idx);
            end
            this.neuronpos = round(r/sys.h(1)) + 1;
        end
        
        function initDynamicAmpPS(this, opts)
            if opts.DynamicPropagationSpeed || opts.DynamicAmplitudes
                base = this.System.Model.DataDir;
                datafile = fullfile(base,'amplitudes.mat');
                s = load(datafile);
                this.amp_xiscale = s.ximax;
                this.amp_kexp = s.amp;
                datafile = fullfile(base,'propagationspeed.mat');
                s = load(datafile);
                this.ps_xiscale = s.ximax;
                this.ps_kexp = s.ps;
                s = load(models.motoneuron.Model.FILE_UPPERLIMITPOLY);
                this.upperlimit_poly = s.upperlimit_poly;
            end
        end
        
    end
    
    
    %% componentwise evaluation
    methods(Access=protected)
        function fx = evaluateComponents(this, pts, ends, idx, self, x, t, mu)
            % This is the template method that actually evaluates the components at given points
            % and values.
            %
            % @attention This method must be able to handle vector-arguments
            % for `\vx,t,\vmu`!
            %
            % Parameters:
            % pts: The components of `\vf` for which derivatives are required @type rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1)+1{:}ends(i)));` @type rowvec<integer>
            % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % x: A matrix `\vX` with the state space locations `\vx_i` in
            % its columns. In rows end(i-1)+1:end(i) it contains the states
            % relevant to evaluate the i-th component of �\vf�.
            % States occur multiply in �\vx� if different components
            % of �\vf� depend on these states.
            % @type matrix<double>
            % t: The corresponding times `t_i` for each state `\vx_i` @type rowvec<double>
            % mu: The corresponding parameters `\mu_i` for each state `\vx_i`, as column matrix
            % @type matrix<double>
            %
            % Return values:
            % fx: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
            % many columns as `\vX` had.
            error('todo');
            fx = zeros(length(pts), length(t));
            
            [j, i, s] = find(this.B(pts,:)');
            idx=zeros(3,length(j));
            idx(3,:) = ceil(j'/this.dims(1)/this.dims(2));
            idx(2,:) = ceil((j'-(idx(3,:)-1)*this.dims(1)*this.dims(2))/this.dims(1));
            idx(1,:) = j'-this.dims(1)*this.dims(2)*(idx(3,:)-1)-this.dims(1)*(idx(2,:)-1);
            neurpos = this.neuronpos;
            MUdistr = this.MUGeoTypeIdx;
            MUActivation = this.computeMUActivation(t);
            ii = MUdistr(idx(2,:)+(idx(3,:)-1)*this.dims(2));
            %jj = neurpos(jdx,kdx):this.dims(1)-1+neurpos(jdx,kdx);
            jj = neurpos(idx(2,:)+(idx(3,:)-1)*this.dims(2))+idx(1,:)-1;
            
            % stretched variant of B
            B2i = i;
            B2j = 1:length(s);
            B2s = s;
            B2 = sparse(B2i, B2j, B2s, length(pts), length(s));
            nmu = length(this.MUTypes);
            for tdx = 1:length(t)
                %act = MUActivation(MUdistr(idx(2,:)+(idx(3,:)-1)*this.dims(2)), neurpos(jdx,kdx):this.dims(1)-1+neurpos(jdx,kdx), t(tdx));       % pseudocode
                act = MUActivation(ii + (jj-1)*nmu + (tdx-1)*nmu*(2*this.dims(1)-1));
                % the used way of getting i, js and s and constructing act
                % ensures that the entries are in the right order
                fx(:,tdx) = B2*act';
            end
        end
        
        function dfx = evaluateComponentPartialDerivatives(this, pts, ends, idx, deriv, self, x, t, mu, dfxsel)%#ok
            % Computes specified partial derivatives of `f` of the components given by pts and
            % the selected partial derivatives by dfxsel.
            %
            % Parameters:
            % pts: The components of `f` for which derivatives are required @type
            % rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % deriv: The indices within `\vx` that derivatives are required for.
            % @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % x: The state space location `\vx` @type colvec<double>
            % t: The corresponding times `t` for the state `\vx` @type double
            % mu: The corresponding parameter `\mu` for the state `\vx` @type colvec<double>
            % dfxsel: A derivative selection matrix. Contains the mapping for each row of x to
            % the output points pts. As deriv might contain less than 'size(x,1)' values, use
            % 'dfxsel(:,deriv)' to select the mapping for the actually computed derivatives.
            %
            % Return values:
            % dfx: A column vector with 'numel(deriv)' rows containing the derivatives at all
            % specified pts i with respect to the coordinates given by 'idx(ends(i-1):ends(i))'
            %
            % See also: setPointSet
            
        end
        
        function fx = evaluateComponentsMulti(this, pts, ends, idx, self, x, t, mu)
            % @todo improve performance!
            fx = zeros(length(pts),length(t));
            mu = mu(1:this.System.ParamCount,:);
            if all(all(repmat(mu(:,1),1,size(mu,2)-1) == mu(:,2:end)))
                this.prepareSimulation(mu(:,1));
                for k=1:length(t)
                    fx(:,k) = this.evaluateComponents(pts, ends, idx, self, [], t(k), mu(:,k));
                end
            else
                for k=1:length(t)
                    this.prepareSimulation(mu(:,k));
                    fx(:,k) = this.evaluateComponents(pts, ends, idx, self, [], t(k), mu(:,k));
                end
            end
        end
    end
end
