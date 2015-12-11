classdef Dynamics < dscomponents.ACompEvalCoreFun
    % This class implements the nonlinear continuum mechanics as described
    % in @cite Heidlauf2013 .
    
    properties
        % Cross-fibre markert part
        b1cf = 5.316372204148964; % [MPa] = [N/mm²]
        d1cf = 0.014991843974911; % [-]
        
        % The activation of the muscle at time t
        %
        % @type function_handle @default @(t)0
        alpha = @(t)0; % [-]
        
        %% Unassembled stuff
        ComputeUnassembled = false;
        % Sigma assembly matrix
        Sigma;
        % The indices of any dirichlet value in the unassembled vector duvw
        idx_uv_bc_glob_unass;
        fDim_unass;
        num_dv_unass;
        idx_v_elems_unass;
        
        %% Force length function fields
        ForceLengthFun;
        ForceLengthFunDeriv;
    end
    
    properties(SetAccess=private)
        APExp;
        AnisoPassiveMuscle;
        AnisoPassiveMuscleDeriv;
        AnisoPassiveTendon;
        AnisoPassiveTendonDeriv;
    end
    
    properties(SetAccess=protected)
        % Helper variable for fullmodels.muscle.model
        lambda_dot_pos;
        lambda_dot;
        
        nfibres;
    end
    
    properties(Transient, SetAccess=private)
        % Prepared arguments for APExpansion
        %         muprep;
        LastBCResiduals;
    end
    
    properties(Transient)
        % Helper value for QuickReleaseTests (or others) that use a
        % function handle with certain alpha ramp time for different
        % simulations. Used in getSimCacheExtra to uniquely identify a
        % simulation in the cache
        RampTime;
    end
    
    properties(Transient, Access=private)
        % Cached quantity from this.fsys.UseDirectMassInversion for
        % faster evaluation of dynamics.
        usemassinv;
        
        % Cached value for cross fibre computations (speed)
        crossfibres = false;
    end
    
    properties(Access=private)
        % Reference to the full system
        fsys;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@dscomponents.ACompEvalCoreFun(sys);
            
            %% Load AP expansion
            %             d = fileparts(which('models.muscle.Dynamics'));
            %             s = load(fullfile(d,'AP'));
            %             s.kexp.Ma = s.kexp.Ma(1,:);
            %             this.APExp = s.kexp;
        end
        
        function configUpdated(this)
            sys = this.fsys;
            mc = sys.Model.Config;
%             if ~isempty(mc.Pool)
%                 this.nfibres = size(mc.FibreTypeWeights,2);
%             end
            if ~isempty(mc)
                this.xDim = sys.NumTotalDofs;
                this.fDim = sys.NumDerivativeDofs;
                this.JSparsityPattern = this.computeSparsityPattern;
                
                %% Sigma assembly matrix
                this.precomputeUnassembledData;
            end
        end
        
        function prepareSimulation(this, mu)
            prepareSimulation@dscomponents.ACompEvalCoreFun(this, mu);
            
            sys = this.fsys;
            mc = sys.Model.Config;
            if ~isempty(mc.Pool)
                mc.Pool.prepare(mu(4),sys.Model.T,sys.Model.dt);
                % Ensure we have a normalized force generated - otherwise
                % most simulations break at 15-fold force
                mc.Pool.NormalizedForce = true;
            end
            % Returns an all zero function if mu(2) is less or equal to zero!
            this.alpha = mc.getAlphaRamp(mu(2));
            
            % Prepare force-length fun
            mc.setForceLengthFun(this);
            
            % Muscle anisotropic passive law
            mlfg = mc.getAnisoMuscleLaw(mu(5),mu(6));
            [this.AnisoPassiveMuscle, this.AnisoPassiveMuscleDeriv] = mlfg.getFunction;
            
            % Tendon anisotropic passive law
            mlfg = mc.getAnisoTendonLaw(mu(7),mu(8));
            [this.AnisoPassiveTendon, this.AnisoPassiveTendonDeriv] = mlfg.getFunction;

            % Get the law function handles that also take b,d as arguments.
            % Needed due to possibly inhomogeneous material.
            %[~,~,this.MarkertLawFun,this.MarkertLawFunDeriv] = mlfg.getFunction;
            
            % Cache stuff
            this.usemassinv = sys.UseDirectMassInversion;
            this.crossfibres = sys.HasFibres && sys.UseCrossFibreStiffness;
        end
        
        function fx = evaluateCoreFun(this, x, t)
            error('Do not call directly; have custom evaluate method.');
        end
        
        function fx = evaluateComponentSet(this, nr, x, t)
            % Manual override of the ACompEvalCoreFun.evaluateComponentSet
            % method as we're doing the "inefficient" approach of computing the
            % full solution for now
            fx = this.evaluate(x,t);
            if ~isempty(this.W)
                fx = this.W(this.PointSets{nr},:)*fx;
            else
                fx = fx(this.PointSets{nr},:);
            end
        end
        
        function fx = evaluateComponentSetMulti(this, nr, x, t, mu)
            % Manual override of the ACompEvalCoreFun.evaluateComponentSetMulti
            % method as we're doing the "inefficient" approach of computing the
            % full solution for now
            fx = this.evaluateMulti(x,t,mu);
            if ~isempty(this.W)
                fx = this.W(this.PointSets{nr},:)*fx;
            else
                fx = fx(this.PointSets{nr},:);
            end
        end
        
        function dfx = evaluateComponentSetGradientsAt(this, nr, x, t)
            % Manual override of the ACompEvalCoreFun.evaluateComponentSetGradientsAt
            % method as we're doing the "inefficient" approach of computing the
            % full solution for now
            J = this.getStateJacobian(x,t);
            if ~isempty(this.W)
                dfx = this.W(this.PointSets{nr},:)*J;
            else
                dfx = J(this.PointSets{nr},:);
            end
        end
        
        function res = test_Jacobian(this, y, t, mu)
            % Overrides the random argument jacobian test as restrictions
            % on the possible x values (detF = 1) hold.
            %
            % Currently the tests using viscosity are commented out as we
            % assume linear damping, which is extracted as extra `A(t,\mu)`
            % part in the models' system
            
            if nargin < 4
                mu = this.System.Model.DefaultMu;
                if nargin < 3
                    t = 1000;
                    if nargin < 2
                        y = this.System.getX0(mu);
                    end
                end
            end
            
            % Use nonzero t to have an effect
            res = test_Jacobian@dscomponents.ACoreFun(this, y, t, mu);
        end
        
        function copy = clone(this)
            % Create new instance
            copy = models.muscle.Dynamics(this.System);
            
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            
            % Copy local properties
            % No local properties are to be copied here, as so far everything is done in the
            % constructor.
        end
        
        function setSystem(this, sys)
            setSystem@dscomponents.ACoreFun(this, sys);
            if isa(sys,'models.ReducedSecondOrderSystem')
                this.fsys = sys.Model.FullModel.System;
            else
                this.fsys = sys;
            end
        end
        
        function copy = project(this, V, W)
             copy = project@dscomponents.ACompEvalCoreFun(this, V, W);
        end
        
        function set.ComputeUnassembled(this, value)
            if value && ~isempty(this.fDim_unass)%#ok
                this.fDim = this.fDim_unass;%#ok
            elseif ~value && ~isempty(this.System)
                this.fDim = this.System.NumDerivativeDofs;
            end
            this.ComputeUnassembled = value;
        end
        
        % Declare public
        [SPK, SPalpha, SPLamDot] = computeSparsityPattern(this);
    end
    
    methods(Access=protected)
        function fx = evaluateComponents(this, pts, ends, ~, ~, x, ~)
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
        end
        
        function dfx = evaluateComponentPartialDerivatives(this, pts, ends, idx, deriv, self, x, t, dfxsel)%#ok
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
    end
    
    methods(Access=private)
        function precomputeUnassembledData(this)
            sys = this.System;
            mc = sys.Model.Config;
            geo = mc.FEM.Geometry;
            this.num_dv_unass = 3*geo.NumElements * geo.DofsPerElement;
            
            % Velocity part: x,y,z velocities
            [i, ~] = find(mc.FEM.Sigma);
            I = [3*(i'-1)+1; 3*(i'-1)+2; 3*(i'-1)+3];
            
            if numel(I) ~= this.num_dv_unass
                error('Size mismatch. Somethings wrong with FEM Sigma');
            end
            
            n = numel(I);
            S = sparse(I,1:n,ones(n,1),3*geo.NumNodes,n);
            
            % Take out nodes with dirichlet BC on output side
            S(sys.idx_v_bc_local,:) = [];
            % Find corresponding unassembled dofs that would be ignored
            % (due to dirichlet velocity values, pressure dirichlet not
            % implemented)
            bc_unass = find(sum(S,1) == 0);
            % Remove them, too. The unassembled evaluation also removes the
            % corresponding entries of the unassembled vector.
            S(:,bc_unass) = [];
            this.idx_uv_bc_glob_unass = bc_unass; %[sys.idx_u_bc_glob' bc_unass];
            this.Sigma = S;
            this.fDim_unass = size(S,2);
            
            % Create boolean array that determines which unassembled dofs
            % belong to which element
            % dK dofs
            hlp = repmat(1:geo.NumElements, 3*geo.DofsPerElement,1);
            hlp = hlp(:);
            hlp(bc_unass) = [];
            ass = false(geo.NumElements,length(hlp));
            for k = 1:geo.NumElements
                ass(k,:) = hlp == k;
            end
            this.idx_v_elems_unass = ass;
        end
    end
end
