classdef RotationFun < dscomponents.ACoreFun
    %ROTATIONDYNSYS Synthetic 2D dynamical system with rotation
    %   Also implements the ACoreFun interface as the target function is
    %   quite simple.
        
    methods
        function this = RotationFun(sys)
            this = this@dscomponents.ACoreFun(sys);
            this.TimeDependent = false;
            this.xDim = 2;
            this.fDim = 2;
        end
        
        function fx = evaluateCoreFun(this, x, t)
            % Implements ACoreFun.evaluate
            a = t*this.mu(1);
            b = a+this.mu(2);
            %b = a+mu(2)+sin(t)/2;
            A = [cos(a) -sin(b); 
                 sin(a) cos(b)];
            fx = A*x;
        end
    end
    
end

