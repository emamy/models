classdef Dynamics < models.muscle.Dynamics;
% Dynamics: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-16
%
% @new{0,7,dw,2014-09-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    %% Properties for motoneuron - sarcomere linking
    properties
        FrequencyDetector        
    end
    
    properties(Dependent)
        UseFrequencyDetector;
    end

    properties(Access=private)
       sarco_force_pos;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@models.muscle.Dynamics(sys);
            
        end
        
        function configUpdated(this)
            sys = this.System;
            mc = sys.Model.Config;
            this.nfibres = length(mc.FibreTypes);
            
            this.lambda_dot_pos = double.empty(2,0);
            this.FrequencyDetector = models.fullmuscle.FrequencyDetector(this.nfibres);
            if sys.HasSpindle
                this.lambda_dot_pos = mc.SpindlePositions;
                this.lambda_dot = zeros(1,this.nfibres);
            end
            
            configUpdated@models.muscle.Dynamics(this);
            
            % fDim and xDim are from models.muscle.Dynamics, so add moto+sarco
            this.fDim = this.fDim + (6+56)*this.nfibres;
            this.xDim = this.xDim + (6+56)*this.nfibres;
            if sys.HasSpindle
                this.fDim = this.fDim + 9*this.nfibres;
                this.xDim = this.xDim + 9*this.nfibres;
            end
            % The forces of the sarcomere in the global vector are in the
            % FirstOrderDofs segment at the positions given by FO
            %
            % Needed in evaluate and getStateJacobian
            this.sarco_force_pos = sys.EndSecondOrderDofs + sys.FO.sarco_output_idx;
        end
        
%         function prepareSimulation(this, mu)
%             prepareSimulation@models.muscle.Dynamics(this, mu);
%         end
        
        function dy = evaluate(this, y, t)
            %% Mechanics
            sys = this.System;
            % Use uvp as argument and also pass in s (=sarco forces)
            force = max(0,y(this.sarco_force_pos)-sys.FO.sarco_mech_signal_offset);
            dy = evaluate@models.muscle.Dynamics(this, y, t, force);
        end
        
        function [SP, SPalpha, SPLamDot] = computeSparsityPattern(this)
            [SP, SPalpha, SPLamDot] = computeSparsityPattern@models.muscle.Dynamics(this);
            % Sarco -> Mechanics
            sys = this.System;
            fo = sys.FO;
            off_sarco = sys.EndSecondOrderDofs + fo.num_motoneuron_dof;
            SP(1:sys.NumDerivativeDofs, off_sarco+(1:fo.num_sarco_dof)) = SPalpha;
        end
        
        function [J, JLamDot] = getStateJacobianImpl(this, y, t)
%             J = this.getStateJacobianFD(y,t);
%             return;
            sys = this.System;
            
            %% Mechanics
            % Use uvp as argument and also pass in s (=sarco forces)
            force = max(0,y(this.sarco_force_pos)-sys.FO.sarco_mech_signal_offset);
            [J, Jalpha, JLamDot] = getStateJacobianImpl@models.muscle.Dynamics(this, y, t, force);
            
            %% Sarcomere -> Mechanics coupling
            % The JS matrix is generated during the computation of the
            % mechanics jacobian, as the element/gauss loop is computed
            % there anyways. its not 100% clean OOP, but is faster.
            fo = sys.FO;
            off_sarco = sys.EndSecondOrderDofs + fo.num_motoneuron_dof;
            J(sys.NumStateDofs + (1:sys.NumDerivativeDofs), off_sarco+(1:fo.num_sarco_dof)) = Jalpha;
        end
        
        function res = test_Jacobian(this, y, t, mu)
            % Overrides the random argument jacobian test as restrictions
            % on the possible x values (detF = 1) hold.
            %
            % Currently the tests using viscosity are commented out as we
            % assume linear damping, which is extracted as extra `A(t,\mu)`
            % part in the models' system
            sys = this.System;
            if nargin < 4
                mu = sys.Model.DefaultMu;
                if nargin < 3
                    t = 1000;
                    if nargin < 2
                        y = this.System.x0.evaluate(mu);
                    end
                end
            end
            
            % Also test correct computation of JLamDot
            d = size(y,1);
            dx = ones(d,1)*sqrt(eps(class(y))).*max(abs(y),1);
            sys.prepareSimulation(mu, sys.Model.DefaultInput);
            this.evaluate(y,t);
            ldotbase = this.lambda_dot;
            uv = sys.NumStateDofs+sys.NumDerivativeDofs;
            LAM = repmat(ldotbase',1,uv);
            dlam = zeros(this.nfibres,uv);
            for k = 1:uv
                h = zeros(d,1);
                h(k) = dx(k);
                this.evaluate(y+h,t);
                dlam(:,k) = this.lambda_dot;
            end
            JLamFD = (dlam - LAM)./dx(1:uv)';
            [~, JLamDot] = this.getStateJacobian(y,t);
            difn = norm(JLamFD-JLamDot);
            res = difn < 1e-13;
            
            freq = ones(1,this.nfibres)*30;
            dx = ones(this.nfibres,1)*sqrt(eps(class(ldotbase))).*max(abs(ldotbase),1);
            sp = sys.Spindle;
            for k = 1:this.nfibres
                spindle_pos = sys.off_spindle + 9*(k-1) + (1:9);
                fx = sp.dydt(y(spindle_pos),t,freq(k),ldotbase(k),0);
                fxh_dldot = sp.dydt(y(spindle_pos),t,freq(k),ldotbase(k)+dx(k),0);
                fxh_dmoto = sp.dydt(y(spindle_pos),t,freq(k)+dx(k),ldotbase(k),0);
                Jspin_Ldot_FD = (fxh_dldot-fx) / dx(k);
                Jspin_moto_FD = (fxh_dmoto-fx) / dx(k);
                [~, Jspin_dLdot, Jspin_dmoto] = sp.Jdydt(y(spindle_pos), t, freq(k), ldotbase(k), 0);
                difn = norm(Jspin_Ldot_FD - Jspin_dLdot');
                res = res && difn < 1e-7;
                difn = norm(Jspin_moto_FD - Jspin_dmoto');
                res = res && difn < 1e-7;
            end
            
            res = res & test_Jacobian@models.muscle.Dynamics(this, y, t, mu);
        end
    end
    
    methods(Static)
        function res = test_Dynamics
            m = models.fullmuscle.Model(fullmuscle.CPull(1));
            f = m.System.f;
            res = true;
            
            f.UseFrequencyDetector = false;
            [t, y] = m.simulate;
            res = res & f.test_Jacobian;
            res = res & f.test_Jacobian(y(:,end),t(end));
            
            f.UseFrequencyDetector = true;
            [t, yfd] = m.simulate;
            res = res & f.test_Jacobian;
            res = res & f.test_Jacobian(yfd(:,end),t(end));
        end
    end
end