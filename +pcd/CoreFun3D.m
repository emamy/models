classdef CoreFun3D < models.pcd.BaseCoreFun
    % The core nonlinear function of the PCD model in 3D.
    %
    % @author Daniel Wirtz @date 2011-10-05
    %
    % @change{0,5,dw,2011-11-02} Augmenting the mu parameters by the base system's
    % models.pcd.BasePCDSystem.ReacCoeff vector. This removes the reaction coefficients from
    % the system as true parameters but allows to quickly revert the process if needed.
    %
    % @new{0,5,dw,2011-10-05} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.morepas.org/software/index.html
    % - \c Documentation http://www.morepas.org/software/kermor/index.html
    % - \c License @ref licensing
    
    properties(SetObservable)
        % The number of seconds after which the neumann conditions decease.
        tDecaySecs = 2;
    end
    
    properties(Access=private)
        A;
        
        nodes;
        
        g;
    end
    
    methods
        
        function this = CoreFun3D(dynsys)
            this = this@models.pcd.BaseCoreFun(dynsys);
        end
        
        function copy = clone(this)
            copy = models.pcd.CoreFun3D(this.sys);
            
            % Call superclass method
            copy = clone@models.pcd.BaseCoreFun(this, copy);
            
            % copy reference!
            %copy.sys = this.sys; % already done in constructor
            copy.A = this.A;
            copy.dim = this.dim;
            copy.g = this.g;
            copy.nodes = this.nodes;
        end
        
        function newSysDimension(this)
            % This fcn is called before each simulation is started.
            
            % Create diffusion matrix
            n = prod(this.System.Dims);
            this.nodes = this.System.geo.Points;
            this.fDim = 4*n;
            this.xDim = 4*n;
        end
        
        function fx = evaluateCoreFun(this, x, t)
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.nodes;
            s = this.System;
            n = this.xDim;
            g = s.geo;
            
            % Uncomment if reaction coeffs become real params again
            mu = [s.ReacCoeff; this.mu];
            
            % Extract single functions
            xa = x(1:m);
            xan = xa.^n;
            ya = x(m+1:2*m);
            xi = x(2*m+1:3*m);
            yi = x(3*m+1:end);
            
            %% Compile boundary conditions
            % Determine the points on each side that are affected of
            % boundary conditions (constant over each simulation as mu is
            % fixed then)
            h = s.h;
            a = s.Omega;
            xr = a(1,2)-a(1,1);
            yr = a(2,2)-a(2,1);
            zr = a(3,2)-a(3,1);
            rb = zeros(m,1);
            ud = (t < this.tDecaySecs) + (t >= this.tDecaySecs)*max(0,2-.5*t);
            
            %% Front & back (only x/y coords relevant)
            [i,j] = ind2sub(g.Dims,g.F); %#ok<*PROP>
            xd = abs((i-1)*h-.5*xr);
            yd = abs((j-1)*h-.5*yr);
            rbp = xd < xr*mu(9)/2 & yd < yr*mu(9)/2;
            rbF = sub2ind(g.Dims, i(rbp), j(rbp), ones(size(i(rbp))));
            rbp = xd < xr*mu(10)/2 & yd < yr*mu(10)/2;
            rbBa = sub2ind(g.Dims, i(rbp), j(rbp), ones(size(i(rbp)))*g.Dims(3));
            rb(rbF) = (xi(rbF)*mu(15)*ud)/h;
            rb(rbBa) = rb(rbBa) + (xi(rbBa)*mu(16)*ud)/h;
            
            %% Left & Right (only x/z coords relevant)
            [i,~,k] = ind2sub(g.Dims,g.L);
            xd = abs((i-1)*h-.5*xr);
            zd = abs((k-1)*h-.5*zr);
            rbp = xd < xr*mu(11)/2 & zd < zr*mu(11)/2;
            rbL = sub2ind(g.Dims, i(rbp), ones(size(i(rbp))), k(rbp));
            rbp = xd < xr*mu(12)/2 & zd < zr*mu(12)/2;
            rbR = sub2ind(g.Dims, i(rbp), ones(size(i(rbp)))*g.Dims(2), k(rbp));
            rb(rbL) = rb(rbL) + (xi(rbL)*mu(17)*ud)/h;
            rb(rbR) = rb(rbR) + (xi(rbR)*mu(18)*ud)/h;
            
            %% Top & Bottom (only y/z coords relevant)
            [i,j,k] = ind2sub(g.Dims,g.T);
            yd = abs((j-1)*h-.5*yr);
            zd = abs((k-1)*h-.5*zr);
            rbp = yd < yr*mu(13)/2 & zd < zr*mu(13)/2;
            rbT = sub2ind(g.Dims, ones(size(i(rbp))), j(rbp), k(rbp));
            rbp = yd < yr*mu(14)/2 & zd < zr*mu(14)/2;
            rbBo = sub2ind(g.Dims, ones(size(i(rbp)))*g.Dims(1), j(rbp), k(rbp));
            rb(rbT) = rb(rbT) + (xi(rbT)*mu(19)*ud)/h;
            rb(rbBo) = rb(rbBo) + (xi(rbBo)*mu(20)*ud)/h;
            
            %% Actual evaluation
            % Xa
            fx(1:m) = mu(1)*xi.*ya - mu(3)*xa + rb;
            % Ya
            fx(m+1:2*m) = mu(2)*yi.*xan - mu(4)*ya;
            % Xi
            fx(2*m+1:3*m) = -mu(1)*xi.*ya - mu(5)*xi + mu(7) - rb;
            % Yi
            fx(3*m+1:end) = -mu(2)*yi.*xan - mu(6)*yi + mu(8);
        end
    end
    
    methods(Access=protected)
        function fxj = evaluateComponents(this, pts, ends, ~, ~, X, t)
            error('Not yet implemented');
        end
    end
end

