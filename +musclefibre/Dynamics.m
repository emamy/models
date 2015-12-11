classdef Dynamics < models.motorunit.MotorunitBaseDynamics
    % Dynamics: Class for nonlinear dynamics of muscle fibre compound
    %
    % Vectorial evaluation of this function is NOT in the sense that you can pass arbitrary values of `t,y,\mu`
    % to evaluate(), but is made such that you can start multiple solves of the system for multiple initial values
    % at the same time. However, this feature is not yet implemented in KerMor.
    % This is due to the used FrequencyDetectors, which need successive time steps in order to work.
    %
    % @author Daniel Wirtz @date 2012-11-22
    %
    % @new{0,8,dw,2015-09-15} Imported into +models package.
    %
    % @new{0,7,dw,2012-11-22} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties(Access=private)
        % Signal to Frequency Converter (linking neuro to spindle)
        %
        % See also: musclefibres.SignaltoFrequencyConverter
        fd;
        options;
    end
    
    methods
        function this = Dynamics(sys, options)
            this = this@models.motorunit.MotorunitBaseDynamics(sys);
            this.options = options;
            
            if this.options.Spindle
                this.fd = SignaltoFrequencyConverter(sys.Model,50,30,2);
            end
            this.fDim = sys.N*sys.dsa + sys.ds + sys.dm;
            this.xDim = this.fDim;
            this.MSLink.MSLink_MaxFactor = 0.25/sys.dx;
            
            this.initJSparsityPattern;
        end
        
        function delete(this)
            this.fd = [];
        end
        
        function dy = evaluate(this, y, t, mu)
            % Evaluates the nonlinear core function at given time and
            % state. Can evaluate vectorized arguments. Here, each column
            % represents one state.
            % The parts for motoneuron, spindle and sarcomeres are evaluated seperately in
            % respective rate functions and the rates are concatenated in this function. The sarcomere part is reshaped such
            % that a matrix with 58 rows and one column for each single sarcomere is passed to the function
            % for (single) sarcomere dynamics.
            %
            % Parameters:
            % y: @type matrix<double>, each column is one state vector
            % t: current time @type rowvec<double>
            % mu: fibre type parameter @type rowvec<double> 0 = slow twitch fibre, 1 = fast twitch fibre
            
            sys = this.System;
            dm = sys.dm;
            ds = sys.ds;
            dsa = sys.dsa;
            N = sys.N;
            
            % Need only fibre_type parameter in here
            mu = this.mu(1,:);
            
            %% Inner dynamics (with partial nonlinear linking)
            dy = zeros(size(y));
            
            % Neuro dynamics (last 2 arguments are 2 phases: phase_prim, phase_sec)
            if this.options.Spindle
                dy(1:dm,:) = sys.moto.dydt(y(1:dm,:), t, y(dm+10:dm+11,:));
                
                % Spindle dynamics
                dy(dm+1:dm+ds,:) = sys.spindle.dydt(y(dm+1:dm+ds,:), t, mu, y(1:dm,:));
            else
                dy(1:dm,:) = sys.moto.dydt(y(1:dm,:), t);
            end
            
            % Sarcomere dynamics
            sa = sys.sarco;
            sadydt = sa.dydt(reshape(y(dm+ds+1:end,:),dsa,[]), t);
            dy(dm+ds+1:end,:) = reshape(sadydt, dsa*N,[]);
            
            %% Link of motoneuron to sarcomer cell
            if ~sys.HadPeak
                linkpos = sys.MotoSarcoLinkIndex;
                dy(linkpos,:) = dy(linkpos,:) + this.mslinkfun(y(2,:))...
                    .*y(2,:)./sa.SarcoConst(1,:);
            end
        end
        
        function dy = evaluateCoreFun(this, y, t, mu)%#ok
            error('evaluate is overridden directly.');
        end
        
        function J = getStateJacobian(this, y, t)
            % Evaluates the full Jacobian of the nonlinear core function.
            % Vectorized evaluation is not tested yet.
            %
            % Parameters:
            % y: @type matrix<double>, each column is one state vector
            % t: current time @type rowvec<double>
            
            s = this.System;
            dm = s.dm;
            ds = s.ds;
            dsa = s.dsa;
            N = s.N;
            
            if this.options.Spindle
                Jm = s.moto.Jdydt(y(1:dm,:), t, 1, y(dm+10:dm+11,:));
                Jsp = s.spindle.Jdydt(y(dm+1:dm+ds,:), t);
            else
                Jm = s.moto.Jdydt(y(1:dm,:), t, 1);
                Jsp = sparse(0,0);
            end
            sa = s.sarco;
            Jsa = cell(1,N);
            ycolwise = reshape(y(dm+ds+1:end,:),dsa,[]);
            [i,j] = find(sa.JSparsityPattern');
            % New, faster evaluation! Stick all into the Jdydt function
            % and create N sparse matrices from that!
            % Tricky pattern transpose and swap of i,j though
            Jsaall = sa.Jdydt(ycolwise, t);
            for idx = 1:N
%                 Jsa{idx} = sa.Jdydt(ycolwise(:,idx), t);
                Jsa{idx} = sparse(j,i,Jsaall(:,idx),dsa,dsa);
            end
            J = blkdiag(Jm,Jsp,Jsa{:});
            J(s.MotoSarcoLinkIndex,2) = this.mslinkfun(y(2,:))./sa.SarcoConst(1,:);
        end
        
        function copy = clone(this)
            % Create new instance
            copy = Dynamics(this.System);
            
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            
            % Copy local properties
            % No local properties are to be copied here, as so far everything is done in the
            % constructor.
        end
        
        function initJSparsityPattern(this)
            % Initializes the Sparsity pattern of the Jacobian
            %
            % TODO: feedback spindle to neuro
            s = this.System;
            
            % Attention: Using double matrix patterns for blkdiag is way
            % faster than sparse logicals, see blkdiag implementation.
            Jm = double(s.moto.JSparsityPattern);
            Jsp = sparse(0,0);
            if this.options.Spindle
                Jsp = double(s.spindle.JSparsityPattern);
            end
            Jsa = cell(1,s.N);
            for idx = 1:s.N
                Jsa{idx} = double(s.sarco.JSparsityPattern);
            end
            J = blkdiag(Jm,Jsp,Jsa{:});
            J(s.MotoSarcoLinkIndex,2) = true;
            this.JSparsityPattern = logical(J);
        end
    end
    
    methods(Access=protected)
        function fx = evaluateComponents(this, pts, ends, ~, ~, x, ~, mu)
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
            % relevant to evaluate the i-th component of `\vf`.
            % States occur multiply in `\vx` if different components
            % of `\vf` depend on these states.
            % @type matrix<double>
            % t: The corresponding times `t_i` for each state `\vx_i` @type rowvec<double>
            % mu: The corresponding parameters `\mu_i` for each state `\vx_i`, as column matrix
            % @type matrix<double>
            %
            % Return values:
            % fx: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
            % many columns as `\vX` had.
            
            error('todo: update');
            % load Rates
            fn = this.f_neuro;
            fsp = this.f_spindle;
            fs = this.f_sarco;
            
            % Need only fibre_type parameter in here
            mu = mu(1,:);
            
            % compute "constants" that depend on mu, the independent
            % constants are already initialized in initxxxRatesCell
            c_moto = this.getMotoConst(mu);
            c1 = this.SarcoConst_slow;
            c2 = this.SarcoConst_fast;
            c_sarco_dep = (1-mu')*c1 + mu'*c2;            
            dm = this.System.dm; ds = this.System.ds; dsa = this.System.dsa;
            linkpos = dm+ds+(this.System.MotoSarcoLinkIndex-1)*dsa+1; % for link motoneuron -> middle sarcomere
            
            % compute fx
            fx = zeros(length(pts), size(x,2));
            for i = 1:length(pts)
                if i == 1
                    st = 0;
                else
                    st = ends(i-1);
                end
                xidx = st+1:ends(i);
                p = pts(i);
                if p <= dm   % p belongs to motoneuron model
                    % if p so, dass Feedback von spindle dort ankommt, dann
                    % zus`tzliches Argument in f_neuro TODO
                    fx(i,:) = fn{p}(x(xidx,:), c_moto);
                elseif p <= dm + ds  % p belongs to spindle model
                    pp = p - dm;
                    if pp == 1
                        % compute link neuro -> spindle with frequency detector TODO
                        %             signal = y_moto(1,:);
                        %             gamma_dyn=this.fd.updateAndGetCurrentFrequency(t, signal);
                        %             gamma_dyn=0;
                        %             gamma_sta=gamma_dyn;
                        gamma_dyn = 20;
                        fx(i,:) = fsp{pp}(x(xidx,:), gamma_dyn);
                    elseif pp == 2 || pp == 5
                        % compute link neuro -> spindle with frequency detector TODO
                        %             signal = y_moto(1,:);
                        %             gamma_dyn=this.fd.updateAndGetCurrentFrequency(t, signal);
                        %             % gamma_dyn=0;
                        %             gamma_sta=gamma_dyn;
                        gamma_sta = 20;
                        fx(i,:) = fsp{pp}(x(xidx,:), gamma_sta);
                    else
                        fx(i,:) = fsp{pp}(x(xidx,:));
                    end
                else  % p belongs to sarcomere model
                    pp = p - (ceil((p-dm-ds)/dsa)-1)*dsa - dm - ds;  % sarcomere component, e.g. the (dm+ds+1)th component is the first component of the first sarcomere
                    if p == linkpos   % p is membrane voltage of the sarcomere which receives motoneuron input
                        input = x(xidx(1),:);    % motoneuron input, neuron states are the first of the full states => xidx(1)
                        fx(i,:) = fs{pp}(x(xidx(2:end),:)', c_sarco_dep)';
                        fac = this.MSLink_MaxFactor*ones(1,size(input,2));
                        dynfac = input < this.MSLink_MaxFactorSignal;
                        fac(dynfac) = this.getLinkFactor(input(dynfac));
                        fx(i,:) = fx(i,:) + fac.*input./((1-mu)*this.SarcoConst_slow(1) + mu*this.SarcoConst_fast(1));
                    else
                        fx(i,:) = fs{pp}(x(xidx,:)', c_sarco_dep)';
                    end
                    
                    
                end
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
    end
end
