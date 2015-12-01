classdef System < models.BaseFirstOrderSystem
    % System: The global dynamical system used within the MotoModel
    %
    %
    % @author Timm Strecker @date 2014-03-24
    %
    % @change{0,8,dw,2015-09-21} Imported into +models package.
    %
    % @new{0,7,ts,2014-03-24} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(Dependent)
        % the conductivity tensor in intracellular space
        %
        % @type matrix<double>
        sigmai;
        
        % the conductivity tensor in extracellular space
        %
        % @type matrix<double>
        sigmae;
        
        % the conductivity tensor in the body (skin/fat)
        %
        % @type matrix<double>
        sigmao;
    end
    
    properties(Access=private)
        % Internal variables with defaults
        sigmai_ = diag([8.93 0.893 0.893]);
        sigmae_ = diag([6.7 6.7 6.7]);
        sigmao_ = diag([.4 .4 .4]);
    end
    
    properties(SetAccess = private)
        % Size of the muscle. vector with 4 entries representing length,
        % width, muscle height and height of the fat/skin layer in [cm].
        %
        % @default [1;1;1;0]
        Geo;
        
        % Discretization stepwidth in cm.
        %
        % @default [0.1; 0.1; 0.1]
        h;
        
        % Dimension of the discretization.
        %
        % @default [10;10;10]
        dim;
        
        % Dimension of the input u
        %
        % @default 1000
        udim;
        
        % number of nodes in z direction lying within the muscle domain
        %
        % @default 10
        zdim_m;
    end
    
    methods
        function this = System(model, opts)
            % Call superclass constructor
            this = this@models.BaseFirstOrderSystem(model);
            this.Geo = opts.Geo;
            this.dim = opts.Dim;
            
            % Parameter mu
            this.addParam('mean current', 4, 'Range', [0 9], 'Desired', 100); 
            
            this.updateDimensions;
            
            this.f = models.emg.RHS(this);
            this.f.init(opts);
            
            % Get diffusion matrix
            this.assembleA;
            % Output of interest - surface of skin part
            this.assembleC;
        end
        
        function assembleA(this)
            ns = this.NumStateDofs;
            musclegrid = general.geometry.RectGrid3D(this.dim(1),this.dim(2),this.dim(3));
            [ie, je, se] = find(MatUtils.generalizedLaplacemat3D(this.h,musclegrid,this.sigmai_+this.sigmae_));
            io = []; jo = []; so = [];
            ib = []; jb = []; sb = [];  % not used if A is constructed this way
            
            % this variant works best:
            if this.Geo(4) ~= 0  % there is a skin layer
                % strategy: create grids for muscle and skin that are too
                % high. create corresponding generalized Laplacians. remove
                % boundary elements between muscle and skin and concatenate Laplacians.

                idxmax = musclegrid.IndexMatrix(end,end,this.zdim_m);
                idx = (ie<=idxmax);
                %ie = ie(1:mmax); je = je(1:mmax); se = se(1:mmax);
                ie = ie(idx); je = je(idx); se = se(idx);
                skingrid = general.geometry.RectGrid3D(this.dim(1),this.dim(2),this.dim(3)-this.zdim_m+1);
                [io, jo, so] = find(MatUtils.generalizedLaplacemat3D(this.h,skingrid,this.sigmao_));
                idxmin = skingrid.IndexMatrix(1,1,2);
                idx = (io >= idxmin);
                %io = io(smin:end); jo = jo(smin:end); so = so(smin:end);
                io = io(idx); jo = jo(idx); so = so(idx);
                io = io + this.udim - idxmin + 1; jo = jo + this.udim - idxmin + 1;
                % ToDo: consider d/dx sigma (because of product rule) - but not
                % differentiable...
            end
         
            % integral over domain
            ii = (ns+1) * ones(1,ns);
            ji = 1:ns;
            si = prod(this.h)*ones(1,ns);
            
            A = sparse([ie' ib' io' ii], [je' jb' jo' ji], [se' -sb' so' si], ns+1,ns);
            this.A = dscomponents.LinearCoreFun(A);
        end        

        function assembleC(this)
            surfacepos = reshape(1:this.NumTotalDofs,this.dim(1),this.dim(2),[]);
            surfacepos = surfacepos(:,:,end);
            C = speye(this.NumTotalDofs);
            idx = true(1,this.NumTotalDofs);
            idx(surfacepos) = false;
            C(idx,:) = [];
            this.C = dscomponents.LinearOutputConv(C);
        end
        
    end
    
    methods(Access=protected)
        function updateDimensions(this)
            d = this.dim;
            geo = this.Geo;
            this.zdim_m = round(geo(3)/(geo(3)+geo(4))*d(3));
            this.udim = d(1)*d(2)*this.zdim_m;
            this.h = [geo(1:2);geo(3)+geo(4)]./(d-1);
            
            this.NumStateDofs = prod(d);
            updateDimensions@models.BaseFirstOrderSystem(this);
        end
    end
    
    %% setters
    methods
        
        function set.sigmae(this, value)
            if ~all(size(value) == [3 3]) || any(any(value-value')) 
               error('sigmae must be a diagonal 3 x 3 matrix.')
            end
            this.sigmae_ = value;
            this.assembleA;
        end
        
        function value = get.sigmae(this)
            value = this.sigmae_;
        end
        
        function set.sigmai(this, value)
            if ~all(size(value) == [3 3]) || any(any(value-value')) 
               error('sigmae must be a diagonal 3 x 3 matrix.')
            end
            this.sigmai_ = value;
            this.assembleA;
            this.f.assembleB;
        end
        
        function value = get.sigmai(this)
            value = this.sigmai_;
        end
        
        function set.sigmao(this, value)
            if ~all(size(value) == [3 3]) || any(any(value-value')) 
               error('sigmao must be a diagonal 3 x 3 matrix.')
            end
            this.sigmao_ = value;
            this.assembleA;
        end
        
        function value = get.sigmao(this)
            value = this.sigmao_;
        end
        
    end
end

