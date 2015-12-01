classdef Dynamics < dscomponents.ACoreFun
    % Dynamics: Class for nonlinear dynamics of motoneuron
    %
    % @author Daniel Wirtz @date 2012-11-22
    %
    % @new{0,7,dw,2012-11-22} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
   
    properties(SetAccess=private)
        moto;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@dscomponents.ACoreFun(sys);
            this.moto = models.motoneuron.Motoneuron;
            this.fDim = 6;
            this.xDim = 6;
        end
        
        function prepareSimulation(this, mu)
            this.moto.setType(mu(1));
            if this.System.Model.FibreTypeDepMaxMeanCurrent
                mu(2) = this.moto.checkMeanCurrent(mu(2));
            end
            prepareSimulation@dscomponents.ACoreFun(this, mu);
        end
        
        function dy = evaluate(this, y, ~)
            % Parameters:
            % y: @type matrix<double>, each column is one state vector
            % t: current time @type rowvec<double>
            dy = this.moto.dydt(y);
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
            % mu: fibre type parameter, 0 = slow twitch fibre, 1 = fast twitch fibre, @type rowvec<double>
            J = this.moto.Jdydt(y, t, 1);
        end
        
        function copy = clone(this)
            % Create new instance
            copy = models.motoneuron.Dynamics(this.System);
            
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACoreFun(this, copy);
            
            % Copy local properties
            copy.moto = this.moto.clone;
        end
    end
end
