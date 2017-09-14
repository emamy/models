classdef System < models.BaseFirstOrderSystem
% System: The global dynamical system used within the Model
%
% Contains the Dynamics, Linear Diffusion and InputConv components as well as a
% ConstInitialValue.
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
        Dimension;
        A_m;
        C_m;
        Sigma_eff;
        Dx;
        L;
        D;
        initialValue;
    end
    
    methods
        function this = System(model)
            % Call superclass constructor
            this = this@models.BaseFirstOrderSystem(model);
            
            % Set core function
            this.f = models.musclefibre_hh.Dynamics(this);
        end
        
        function newDim(this)
            this.updateDimensions;
        end
    end
    
    methods(Static)
        function [STATES] = initStates()
            %V m h n
            STATES(1) = -75;
            STATES(2) = 0.05;
            STATES(3) = 0.6;
            STATES(4) = 0.325;
        end
        function [KK]= StiffnessMatrix_Vmnh(dim,time_step,dx,D)
            % stiffness matrix (V m h n V m h n)
            n=dim;
            KK = zeros(dim,dim);
  
            % introduce factor for compact notation      
            fact = time_step/dx^2 * D;      
      
            for i=2:n-7
            % Wenn Index = V
            if mod(i, 4) == 1
                KK(i,i-4) = fact;
                KK(i,i)   = -2*fact;
                KK(i,i+4) = fact;    
            end
            end
  
        %   for j = 1:n
        %       if mod(j,4) ~= 1
        %           KK(j,j)  = 1;
        %       end
        %   end
  
         % homogeneous Neumann B.C. - no flow
        % right side
         %KK(n-3,n-3) = -1;
         KK(n-3,n-7) = 1;
        % left side
            %KK(1,1)     = -1;
         KK(1,5)     = 1;  
        end
        
    end
    
    methods(Access=protected)
        
        function updateDimensions(this)
            % um Dimension zu bekommen, etc.
            m = this.Model;
            this.Dimension = m.Dimension;
            this.A_m = m.A_m;
            this.C_m = m.C_m;
            this.Sigma_eff = m.Sigma_eff;
            this.Dx = m.Dx;
            this.L = m.L;
            this.D = this.Sigma_eff/this.A_m/this.C_m;
            % Parameter für x0
            x0 = []';
            for i=1: floor(this.Dimension/4)
                x0= [x0 this.initStates];
            end
            this.initialValue = x0';
            this.x0 = dscomponents.ConstInitialValue(x0');           
            
            % x' = A * x + B * u(t) + f(x, t)
            % vllt. nachher andere Funktion
            this.assembleA;
            %this.A = dscomponents.LinearCoreFun(A);
            
            this.assembleB;
            
            this.assembleC;
            
            this.NumStateDofs = this.Dimension;
            this.f.newDim;            
            updateDimensions@models.BaseFirstOrderSystem(this);
        end
        
        function A = assembleA(this)
            % Computes the linear term A which represents the diffusion of
            % the membrane voltage along the sarcomeres. Diffusion is
            % discretized via finite differences/ central difference.
            A = this.StiffnessMatrix_Vmnh(this.Dimension, this.Model.dt, this.Dx, this.D);
            a = dscomponents.LinearCoreFun(A);
            this.A = a;
           
        
        end  
        
        function B = assembleB(this)
            b = zeros(this.Dimension,floor(this.Dimension/4));
            j=1;
            for i=1:this.Dimension
                if mod(i,4)==1
                    b(i,j)=1;
                    j = j+1;
                end
            end
            b = (1/this.C_m) * b;
            B = dscomponents.LinearInputConv(b);
            this.B = B;
        end
        
        function C = assembleC(this)
            v = floor(this.Dimension / 4);
            c = zeros(v, this.Dimension);
            j=1;
            for i=1:v
                c(i,j) = 1;
                j = j+4;
            end
            C = dscomponents.LinearOutputConv(c);
            this.C = C;
        end
    end
end
