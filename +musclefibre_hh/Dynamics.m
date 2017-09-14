classdef Dynamics < dscomponents.ACompEvalCoreFun
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
    
    properties
        
    end
    
    methods
        function this = Dynamics(sys)
            this = this@dscomponents.ACompEvalCoreFun(sys);
            this.CustomProjection = false;
            this.TimeDependent = false;
        end
        
        function fx = evaluateCoreFun(this, x, ~, ~)
            fx = [];
            for i=1:4:length(x)
                states = [x(i) x(i+1) x(i+2) x(i+3)]';
                fxj = this.computeRates_HODGKIN_HUXLEY(states, this.initConsts);
                fx = [fx ; fxj];
            end
        end
        
        function J = getStateJacobian(this, x, ~, ~)
            J = zeros(length(x));
            constants = this.initConsts;
            for i=1:length(x)
                if mod(i,4) == 1
                    J(i,i:i+3)   = [this.computeRateD_DV_V(x(i:i+3), constants) this.computeRateD_DM_V(x(i:i+3), constants) this.computeRateD_DH_V(x(i:i+3), constants) this.computeRateD_DN_V(x(i:i+3), constants)];
                    J(i+1,i:i+3) = [this.computeRateD_DV_M(x(i:i+3), constants) this.computeRateD_DM_M(x(i:i+3), constants) this.computeRateD_DH_M(x(i:i+3), constants) this.computeRateD_DN_M(x(i:i+3), constants)];
                    J(i+2,i:i+3) = [this.computeRateD_DV_H(x(i:i+3), constants) this.computeRateD_DM_H(x(i:i+3), constants) this.computeRateD_DH_H(x(i:i+3), constants) this.computeRateD_DN_H(x(i:i+3), constants)];
                    J(i+3,i:i+3) = [this.computeRateD_DV_N(x(i:i+3), constants) this.computeRateD_DM_N(x(i:i+3), constants) this.computeRateD_DH_N(x(i:i+3), constants) this.computeRateD_DN_N(x(i:i+3), constants)];
                end
            end
        end
        
        function newDim(this)
            m = this.System.Model;
            n = m.Dimension;
            this.xDim = n;
            this.fDim = n;
            StateJacobian = this.getStateJacobian(this.System.initialValue, []);
            this.JSparsityPattern = spones(sparse(StateJacobian));
        end
        
        function copy = clone(this)
            copy = clone@dscomponents.ACompEvalCoreFun(this, models.musclefibre_hh.Dynamics(this.System));
        end
    end
    
    methods(Access=protected)
        function fxj = evaluateComponents(this, pts, ends, argidx, self, X, ~, ~)
            % Evaluates the burgers nonlinearity pointwise.
            % WIRKLICH 4?
            % Parameters:
            % pts: The components of `\vf` for which derivatives are required @type rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % argidx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % X: A matrix `\vX` with the state space locations `\vx_i` in its columns @type
            % matrix<double>
            %
            % Return values:
            % fxj: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
            % many columns as `\vX` had.
            %anzahl punkte, anzahl zeitschritte
            % für yi, die benötigt werden
            fxj = zeros(length(pts),size(X,2));
            for idx=1:length(pts)
                constants = this.initConsts;
                pt = pts(idx);
                if idx == 1
                    st = 0;
                else
                    st = ends(idx-1);
                end
                % 4, wenn V. Sonst 2(m h n)
                % funktioniert, da sich Werte wiederholen in X, deshalb
                % sind sie aufeinanderfolgend
                % lokaler Index
                xidx = (st+1):ends(idx);
                x = X(xidx,:);
                globidx = argidx(xidx);
                richtigx = globidx * self(xidx)';
                if mod(pt, 4) ==1
                    temp= this.computeRates_HODGKIN_HUXLEY(x, constants);
                    % V wird für alle t_i zurückgegeben
                    fxj(idx,:) = temp(1,:);
                elseif mod(pt, 4) ==2
                    fxj(idx,:)= this.computeRate_M(x, constants);
                elseif mod(pt, 4) ==3
                    fxj(idx,:)= this.computeRate_H(x, constants);
                else
                    fxj(idx,:)= this.computeRate_N(x, constants);
                end
            end
        end
        
        function fxj = evaluateComponentsMulti(this, varargin)
            fxj = this.evaluateComponents(varargin{:});
        end
        
        function dfx = evaluateComponentPartialDerivatives(this, pts, ends, argidx, deriv, self, x, t, dfxsel)
            % Computes specified partial derivatives of `f` of the components given by pts and
            % the selected partial derivatives by dfxsel.
            %
            % Override in subclasses for optimized performance
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
            % dfxsel: A derivative selection matrix. Contains the mapping for each row of x to
            % the output points pts. As deriv might contain less than 'size(x,1)' values, use
            % 'dfxsel(:,deriv)' to select the mapping for the actually computed derivatives.
            %
            % Return values:
            % dfx: A column vector with 'numel(deriv)' rows containing the derivatives at all
            % specified pts i with respect to the coordinates given by 'idx(ends(i-1):ends(i))'
            %
            % See also: setPointSet
            constants = this.initConsts;
            for i=1:length(deriv)
                for idx=1:length(pts)
                    pt = pts(idx);
                    if idx == 1
                        st = 0;
                    else
                        st = ends(idx-1);
                    end
                    % 4, wenn V. Sonst 2(m h n)
                    % funktioniert, da sich Werte wiederholen in X, deshalb
                    % sind sie aufeinanderfolgend
                    % lokaler Index
                    xidx = (st+1):ends(idx);
                    if ismember(deriv(i),xidx)
                        x_eval = x(xidx);
                        derivative = argidx(deriv(i));
                        temp = this.computeRates(x_eval, constants, derivative, pt);
                        dfx(i, :) = temp(1, :);
                    end
                end
            end
        end
        
        function RATES = computeRates(this, x_eval, constants, deriv, pt)
            %x_eval muss je nach Antwort modifiziert werden
            if mod(pt, 4) == 1
                x = x_eval;
                if mod(deriv, 4) == 1
                    RATES = this.computeRateD_DV_V(x, constants);
                elseif mod(deriv, 4) == 2
                    RATES = this.computeRateD_DM_V(x, constants);
                elseif mod(deriv, 4) == 3
                    RATES = this.computeRateD_DH_V(x, constants);
                else
                    RATES = this.computeRateD_DN_V(x, constants);
                end
            elseif mod(pt, 4) == 2
                x = [x_eval(1) x_eval(2) 0 0]';
                if mod(deriv, 4) == 1
                    RATES = this.computeRateD_DV_M(x, constants);
                elseif mod(deriv, 4) == 2
                    RATES = this.computeRateD_DM_M(x, constants);
                elseif mod(deriv, 4) == 3
                    RATES = this.computeRateD_DH_M(x, constants);
                else
                    RATES = this.computeRateD_DN_M(x, constants);
                end
            elseif mod(pt, 4) == 3
                x = [x_eval(1) 0 x_eval(2) 0]';
                if mod(deriv, 4) == 1
                    RATES = this.computeRateD_DV_H(x, constants);
                elseif mod(deriv, 4) == 2
                    RATES = this.computeRateD_DM_H(x, constants);
                elseif mod(deriv, 4) == 3
                    RATES = this.computeRateD_DH_H(x, constants);
                else
                    RATES = this.computeRateD_DN_H(x, constants);
                end
            else
                x = [x_eval(1) 0 0 x_eval(2)]';
                if mod(deriv, 4) == 1
                    RATES = this.computeRateD_DV_N(x, constants);
                elseif mod(deriv, 4) == 2
                    RATES = this.computeRateD_DM_N(x, constants);
                elseif mod(deriv, 4) == 3
                    RATES = this.computeRateD_DH_N(x, constants);
                else
                    RATES = this.computeRateD_DN_N(x, constants);
                end
            end
        end
    end
    
    methods(Static)
        function [CONSTANTS] = initConsts()
            CONSTANTS = [];
            %E_r
            CONSTANTS(:,1) = -75.000000;
            %C_m
            CONSTANTS(:,2) = 1.000000;
            %g_Na
            CONSTANTS(:,3) = 120.000000;
            %g_K
            CONSTANTS(:,4) = 36.000000;
            %g_L
            CONSTANTS(:,5) = 0.3000000;
            %E_Na
            CONSTANTS(:,6) = CONSTANTS(:,1)+115.000000;
            %E_K
            CONSTANTS(:,7) = CONSTANTS(:,1) - 12.000000;
            %E_L
            CONSTANTS(:,8) = CONSTANTS(:,1)+10.613000000;
        end
        
        function [RATES] = computeRates_HODGKIN_HUXLEY(STATES, CONSTANTS)
            statesSize = size(STATES);
            statesColumnCount = statesSize(2);
            if ( statesColumnCount == 1)
                STATES = STATES';
                ALGEBRAIC = zeros(1, 10);
            else
                % davor  t1 t2 t3 t4
                %      V
                %      m
                %      h
                %      n
                STATES = STATES';
                statesSize = size(STATES);
                statesColumnCount = statesSize(2);
                statesRowCount = statesSize(1);
                ALGEBRAIC = zeros(statesRowCount, 10);
                RATES = zeros(statesRowCount, statesColumnCount);
            end
            %alpha_m
            ALGEBRAIC(:,2) = (  - 0.1000000.*(STATES(:,1)+50.000000))./((exp(( - (STATES(:,1)+50.000000)./10.000000))) - 1.000000);
            %beta_m
            ALGEBRAIC(:,6) =  4.000000.*(exp(( - (STATES(:,1)+75.000000)./18.000000)));
            %d/dt m
            RATES(:,2) =  ALGEBRAIC(:,2).*(1.000000 - STATES(:,2)) -  ALGEBRAIC(:,6).*STATES(:,2);
            %alpha_h
            ALGEBRAIC(:,3) =  0.07000000.*(exp(( - (STATES(:,1)+75.000000)./20.000000)));
            %beta_h
            ALGEBRAIC(:,7) = 1.000000./((exp(( - (STATES(:,1)+45.000000)./10.000000)))+1.000000);
            %d/dt h
            RATES(:,3) =  ALGEBRAIC(:,3).*(1.000000 - STATES(:,3)) -  ALGEBRAIC(:,7).*STATES(:,3);
            %alpha_n
            ALGEBRAIC(:,4) = (  - 0.01000000.*(STATES(:,1)+65.000000))./((exp(( - (STATES(:,1)+65.000000)./10.000000))) - 1.000000);
            %beta_n
            ALGEBRAIC(:,8) =  0.125000000.*(exp(((STATES(:,1)+75.000000)./80.000000)));
            %d/dt n
            RATES(:,4) =  ALGEBRAIC(:,4).*(1.00000 - STATES(:,4)) -  ALGEBRAIC(:,8).*STATES(:,4);
            %I_Na
            ALGEBRAIC(:,5) =  CONSTANTS(:,3).*(STATES(:,2) .^ 3.00000).*STATES(:,3).*(STATES(:,1) - CONSTANTS(:,6));
            %I_K
            ALGEBRAIC(:,9) =  CONSTANTS(:,4).*(STATES(:,4) .^ 4.00000).*(STATES(:,1) - CONSTANTS(:,7));
            %I_L
            ALGEBRAIC(:,10) =  CONSTANTS(:,5).*(STATES(:,1) - CONSTANTS(:,8));
            
            %d/dt V
            RATES(:,1) =  - 1.000000 .* ( ALGEBRAIC(:,5)+ALGEBRAIC(:,9)+ALGEBRAIC(:,10))./CONSTANTS(:,2);
            %wieder  t1 t2 t3 t4
            %       V
            %       m
            %       h
            %       n
            RATES = RATES';
        end
        
        function [RATES] = computeRate_M(STATES, ~)
            statesSize = size(STATES);
            statesColumnCount = statesSize(2);
            if ( statesColumnCount == 1)
                STATES = STATES';
                % nur zwei ALGEBRAIC
                ALGEBRAIC = zeros(1, 2);
            else
                % davor  t1 t2 t3 t4
                %      V
                %      m
                %      h
                %      n
                STATES = STATES';
                statesSize = size(STATES);
                statesRowCount = statesSize(1);
                % nur zwei ALGEBRAIC
                ALGEBRAIC = zeros(statesRowCount, 2);
                % nur m wird zurückgegeben, nicht V und m
                RATES = zeros(statesRowCount, 1);
            end
            
            %alpha_m
            ALGEBRAIC(:,1) = (  - 0.1000000.*(STATES(:,1)+50.000000))./((exp(( - (STATES(:,1)+50.000000)./10.000000))) - 1.000000);
            %beta_m
            ALGEBRAIC(:,2) =  4.000000.*(exp(( - (STATES(:,1)+75.000000)./18.000000)));
            %d/dt m
            RATES(:,1) =  ALGEBRAIC(:,1).*(1.000000 - STATES(:,2)) -  ALGEBRAIC(:,2).*STATES(:,2);
            %wieder  t1 t2 t3 t4
            %       V
            %       m
            %       h
            %       n
            RATES = RATES';
        end
        
        function [RATES] = computeRate_H(STATES, ~)
            statesSize = size(STATES);
            statesColumnCount = statesSize(2);
            if ( statesColumnCount == 1)
                STATES = STATES';
                ALGEBRAIC = zeros(1, 2);
            else
                % davor  t1 t2 t3 t4
                %      V
                %      m
                %      h
                %      n
                STATES = STATES';
                statesSize = size(STATES);
                statesRowCount = statesSize(1);
                ALGEBRAIC = zeros(statesRowCount, 2);
                % nur h wird zurückgegeben, nicht V und h
                RATES = zeros(statesRowCount, 1);
            end
            
            %alpha_h
            ALGEBRAIC(:,1) =  0.07000000.*(exp(( - (STATES(:,1)+75.000000)./20.000000)));
            %beta_h
            ALGEBRAIC(:,2) = 1.000000./((exp(( - (STATES(:,1)+45.000000)./10.000000)))+1.000000);
            %d/dt h
            %Rates(:,1), weil nur ein Wert zurückgegeben wird; States(:,2),
            %weil nur 2 Werte übergeben wurden.
            RATES(:,1) =  ALGEBRAIC(:,1).*(1.000000 - STATES(:,2)) -  ALGEBRAIC(:,2).*STATES(:,2);
            %wieder  t1 t2 t3 t4
            %       V
            %       m
            %       h
            %       n
            RATES = RATES';
        end
        
        function [RATES] = computeRate_N(STATES, ~)
            statesSize = size(STATES);
            statesColumnCount = statesSize(2);
            if ( statesColumnCount == 1)
                STATES = STATES';
                ALGEBRAIC = zeros(1, 2);
            else
                % davor  t1 t2 t3 t4
                %      V
                %      m
                %      h
                %      n
                STATES = STATES';
                statesSize = size(STATES);
                statesRowCount = statesSize(1);
                ALGEBRAIC = zeros(statesRowCount, 2);
                % nur n wird zurückgegeben, nicht V und n
                RATES = zeros(statesRowCount, 1);
            end
            
            %alpha_n
            ALGEBRAIC(:,1) = (  - 0.01000000.*(STATES(:,1)+65.000000))./((exp(( - (STATES(:,1)+65.000000)./10.000000))) - 1.000000);
            %beta_n
            ALGEBRAIC(:,2) =  0.125000000.*(exp(((STATES(:,1)+75.000000)./80.000000)));
            %d/dt n
            %Rates(:,1), weil nur ein Wert zurückgegeben wird; States(:,2),
            %weil nur 2 Werte übergeben wurden.
            RATES(:,1) =  ALGEBRAIC(:,1).*(1.000000 - STATES(:,2)) -  ALGEBRAIC(:,2).*STATES(:,2);
            %wieder  t1 t2 t3 t4
            %       V
            %       m
            %       h
            %       n
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DV_V(STATES, CONSTANTS)
            STATES = STATES';
            %d/dv V
            %= -1/C_m * (g_Nam^3h + g_Kn^4+g_L)
            RATES(:,1) =  - ( CONSTANTS(:,3) .* (STATES(:,2).^3.000000) * STATES(:,3) + CONSTANTS(:,4) .* (STATES(:,4).^4.000000) + CONSTANTS(:,5))./CONSTANTS(:,2);
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DM_V(STATES, CONSTANTS)
            STATES = STATES';
            %d/dm V
            %= -1/C_m * (3 * g_Nam^2h(V - E_Na))
            RATES(:,1) =  - (3.000000 .* CONSTANTS(:,3).*(STATES(:,2) .^ 2.000000).*STATES(:,3).*(STATES(:,1) - CONSTANTS(:,6)))./CONSTANTS(:,2);
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DH_V(STATES, CONSTANTS)
            STATES = STATES';
            %d/dh V
            %= -1/C_m * (g_Nam^3(V - E_Na))
            RATES(:,1) =  - (CONSTANTS(:,3).*(STATES(:,2) .^ 3.000000).*(STATES(:,1) - CONSTANTS(:,6)))./CONSTANTS(:,2);
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DN_V(STATES, CONSTANTS)
            STATES = STATES';
            %d/dn V
            %= -1/C_m * (4*g_K*n^3(V - E_K))
            RATES(:,1) =  - (4.000000.*CONSTANTS(:,4).*(STATES(:,4) .^ 3.000000).*(STATES(:,1) - CONSTANTS(:,7)))./CONSTANTS(:,2);
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DV_M(STATES, ~)
            STATES = STATES';
            %d/dv M
            %= d/dv alpha_m * (1-m) - d/dv beta_m*m
            
            %d/dv alpha_m
            %= -0.1 * (-1 + exp( -(V + 50)/10)) - (-0.1 * (V+50)) * exp( -(V +
            %50)/10) * -1/10
            %alles / (-1 + exp( -(V + 50)/10))^2
            ALGEBRAIC(:,1) = ((-0.1000000 .* (-1.000000 + exp(( - (STATES(:,1)+50.000000)./10.000000)))) - ((  - 0.1000000.*(STATES(:,1)+50.000000)) .* exp(( - (STATES(:,1)+50.000000)./10.000000)) .* (-1.000000/10.000000))) ./ (((-1.000000 + exp(( - (STATES(:,1)+50.000000)./10.000000)))).^2.000000);
            
            %d/dv beta_m
            %=4 * exp(-(V+50)/18)*-1/18
            ALGEBRAIC(:,2) =  4.000000.*(exp(( - (STATES(:,1)+75.000000)./18.000000))) .* (-1.000000 / 18.000000);
            
            %d/dt m
            RATES(:,1) =  ALGEBRAIC(:,1).*(1.000000 - STATES(:,2)) -  ALGEBRAIC(:,2).*STATES(:,2);
            
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DM_M(STATES, ~)
            STATES = STATES';
            %d/dm M
            %alpha_m
            ALGEBRAIC(:,1) = (  - 0.1000000.*(STATES(:,1)+50.000000))./((exp(( - (STATES(:,1)+50.000000)./10.000000))) - 1.000000);
            %beta_m
            ALGEBRAIC(:,2) =  4.000000.*(exp(( - (STATES(:,1)+75.000000)./18.000000)));
            %d/dt m
            %= - alpha_m - beta_m
            RATES(:,1) =  - ALGEBRAIC(:,1) -  ALGEBRAIC(:,2);
            
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DH_M(~, ~)
            %d/dh M
            %= 0
            RATES(:,1) =  0.000000;
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DN_M(~, ~)
            %d/dn M
            %= 0
            RATES(:,1) =  0.000000;
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DV_H(STATES, ~)
            STATES = STATES';
            %d/dv H
            %d/dv alpha_h
            %=%0.07 * exp(-(V+75)/20)*(-1/20)
            ALGEBRAIC(:,1) =  0.07000000.*(exp(( - (STATES(:,1)+75.000000)./20.000000))) .* (-(1.000000/20.000000));
            %d/dv beta_h
            %= -(1 + exp(-(V+45)/10))^-2 * exp(-(V +45)/10)
            ALGEBRAIC(:,2) = -(((exp(( - (STATES(:,1)+45.000000)./10.000000)))+1.000000).^-2.000000) .* (exp(( - (STATES(:,1)+45.000000)./10.000000))) .* (-1.000000/10.000000);
            %d/dt h
            %= d/dv alpha_h * (1 - h) - d/dv beta_h * h
            RATES(:,1) =  ALGEBRAIC(:,1).*(1.000000 - STATES(:,3)) -  ALGEBRAIC(:,2).*STATES(:,3);
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DM_H(~, ~)
            %d/dm H
            %= 0
            RATES(:,1) =  0.000000;
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DH_H(STATES, ~)
            STATES = STATES';
            %d/dh H
            %alpha_h
            ALGEBRAIC(:,1) =  0.07000000.*(exp(( - (STATES(:,1)+75.000000)./20.000000)));
            %beta_h
            ALGEBRAIC(:,2) = 1.000000./((exp(( - (STATES(:,1)+45.000000)./10.000000)))+1.000000);
            %d/dt h
            %= -alpha_h - beta_h
            RATES(:,1) =  - ALGEBRAIC(:,1) -  ALGEBRAIC(:,2);
            
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DN_H(~, ~)
            %d/dn H
            %= 0
            RATES(:,1) =  0.000000;
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DV_N(STATES, ~)
            STATES = STATES';
            %d/dv N
            %= d/dv alpha_n * (1-n) - d/dv beta_n * n
            %d/dv alpha_n
            %= -0.01 * (-1+exp(-(V+65)/10)) -
            %(-0.01*(V+65))*exp(-(V+65)/10)*(-1/10)
            %alles/ (-1 + exp(-(V + 65)/10)^2
            ALGEBRAIC(:,1) = (- 0.01000000.* ((exp(( - (STATES(:,1)+65.000000)./10.000000))) - 1.000000) - ((  - 0.01000000.*(STATES(:,1)+65.000000)) .* (exp(( - (STATES(:,1)+65.000000)./10.000000))) .* (-1.000000/10.000000))) ./ ((((exp(( - (STATES(:,1)+65.000000)./10.000000))) - 1.000000)).^2.000000);
            %d/dv beta_n
            %= 0.125 * exp((V+75)/80)*(1/80)
            ALGEBRAIC(:,2) =  0.125000000.*(exp(((STATES(:,1)+75.000000)./80.000000))) .* (1.000000 / 80.000000);
            %d/dt n
            RATES(:,1) =  ALGEBRAIC(:,1).*(1.000000 - STATES(:,4)) -  ALGEBRAIC(:,2).*STATES(:,4);
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DM_N(~, ~)
            %d/dm N
            %= 0
            RATES(:,1) =  0.000000;
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DH_N(~, ~)
            %d/dh N
            %= 0
            RATES(:,1) =  0.000000;
            RATES = RATES';
        end
        
        function [RATES] = computeRateD_DN_N(STATES, ~)
            STATES = STATES';
            %d/dn N
            %alpha_n
            ALGEBRAIC(:,1) =  (  - 0.01000000.*(STATES(:,1)+65.000000))./((exp(( - (STATES(:,1)+65.000000)./10.000000))) - 1.000000);
            %beta_n
            ALGEBRAIC(:,2) = 0.125000000.*(exp(((STATES(:,1)+75.000000)./80.000000)));
            %d/dt n
            %= -alpha_n - beta_n
            RATES(:,1) =  - ALGEBRAIC(:,1) -  ALGEBRAIC(:,2);
            
            RATES = RATES';
        end
    end
end