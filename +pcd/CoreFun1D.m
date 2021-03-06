classdef CoreFun1D < models.pcd.BaseCoreFun
% The core nonlinear function of the PCD model.
%
% @author Daniel Wirtz @date 2010-03-16
%
% @new{0,6,dw,2011-11-26} Computing the JSparsityPattern
%
% @change{0,5,dw,2011-11-02} Augmenting the mu parameters by the base system's
% models.pcd.BasePCDSystem.ReacCoeff vector. This removes the reaction coefficients from
% the system as true parameters but allows to quickly revert the process if needed.
%
% @change{0,3,dw,2011-04-12} Added MultiArgumentEvaluation capabilities to this class.
%
% @new{0,3,dw,2011-03-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(Access=private)
        % The number of nodes
        nodes;
    end
    
    methods
                
        function this = CoreFun1D(dynsys)
            this = this@models.pcd.BaseCoreFun(dynsys);
            this.TimeDependent = false;
        end
        
        function copy = clone(this)
            copy = models.pcd.CoreFun1D(this.System);
            
            % Call superclass method
            copy = clone@models.pcd.BaseCoreFun(this, copy);
            copy.nodes = this.nodes;
        end
        
        function newSysDimension(this)
            n = this.System.Dims(1);
            % Add x_a dependencies
            % 1=x_a, 2=y_a, 3=x_i {, y_i}
            i = [(1:n)'; (1:n)'; (1:n)']; 
            j = (1:3*n)';
            % Add y_a dependencies
            % 1=x_a, 2=y_a {, x_i}, 3=y_i
            i = [i; (n+1:2*n)'; (n+1:2*n)'; (n+1:2*n)'];
            j = [j; (1:2*n)'; (3*n+1:4*n)'];
            % Add x_i dependencies
            % {x_a,} 1=y_a, 2=x_i {, y_i}
            i = [i; (2*n+1:3*n)'; (2*n+1:3*n)']; 
            j = [j; (n+1:3*n)'];
            % Add y_i dependencies
            % 1=x_a, {y_a, x_i}, 2=y_i
            i = [i; (3*n+1:4*n)'; (3*n+1:4*n)']; 
            j = [j; (1:n)'; (3*n+1:4*n)'];
            this.xDim = 4*n;
            this.fDim = 4*n;
            this.JSparsityPattern = sparse(i,j,ones(length(i),1),this.fDim,this.xDim);
            this.nodes = n;
        end
        
        function fx = evaluate(this, x, ~)
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.nodes;
            n = this.System.Model.n;
            
            % Non-Multidimensional case
            mu = [this.System.ReacCoeff; this.mu];

            % Extract single functions
            xa = x(1:m);
            xan = xa.^n;
            ya = x(m+1:2*m);
            xi = x(2*m+1:3*m);
            yi = x(3*m+1:end);

            % Boundary conditions
            rb = zeros(m,1);
            rb(end) = xi(end)*mu(9)/this.System.hs; %max((500-t)/500,0)*

            fx(1:m) = mu(1)*xi.*ya - mu(3)*xa + rb;
            fx(m+1:2*m) = mu(2)*yi.*xan - mu(4)*ya;
            fx(2*m+1:3*m) = -mu(1)*xi.*ya - mu(5)*xi + mu(7) - rb;
            fx(3*m+1:end) = -mu(2)*yi.*xan - mu(6)*yi + mu(8);
        end
        
        function fx = evaluateMulti(this, x, ~, mu)
            % Allocate result vector
            fx = zeros(size(x));
            
            m = this.nodes;
            n = this.System.Model.n;
            
            mu = [repmat(this.System.ReacCoeff,1,size(mu,2)); mu];
            % Extract single functions
            xa = x(1:m,:);
            xan = xa.^n;
            ya = x(m+1:2*m,:);
            xi = x(2*m+1:3*m,:);
            yi = x(3*m+1:end,:);

            % Boundary conditions
            rb = zeros(m,size(xi,2));
            rb(end,:) = (xi(end,:).*mu(9,:))/this.System.hs;

            fx(1:m,:) = bsxfun(@times,xi.*ya,mu(1,:)) - bsxfun(@times,xa,mu(3,:)) + rb;
            fx(m+1:2*m,:) = bsxfun(@times,yi.*xan,mu(2,:)) - bsxfun(@times,ya,mu(4,:));
            fx(2*m+1:3*m,:) = -bsxfun(@times,xi.*ya,mu(1,:)) - bsxfun(@times,xi,mu(5,:)) + bsxfun(@times,ones(size(xi)),mu(7,:)) - rb;
            fx(3*m+1:end,:) = -bsxfun(@times,yi.*xan,mu(2,:)) - bsxfun(@times,yi,mu(6,:)) + bsxfun(@times,ones(size(xi)),mu(8,:));
        end
    end
    
    methods(Access=protected)
        function fxj = evaluateComponents(this, pts, ends, ~, ~, X, ~)
            fxj = this.evaluateComponentsMulti(pts, ends, [], [], X, [], this.mu);
        end
        
        function fxj = evaluateComponentsMulti(this, pts, ends, ~, ~, X, ~, mu)
            % The vector embedding results from the fixed ordering of the full 4*m-vector into
            % the components x_a, y_a, x_i, y_i
            %
            % Parameters:
            % pts: The components of `\vf` for which derivatives are required @type rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % X: A matrix `\vX` with the state space locations `\vx_i` in its columns @type
            % matrix<double>
            % mu: The corresponding parameters `\mu_i` for each state `\vx_i`, as column matrix
            % @type matrix<double>
            %
            % Return values:
            % fxj: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
            % many columns as `\vX` had.
            m = this.nodes;
            s = this.System;
            n = s.Model.n;
            fxj = zeros(length(pts),size(X,2));
            
            mu = [repmat(s.ReacCoeff,1,size(mu,2)); mu];
            for idx=1:length(pts)
                j = pts(idx);
                if idx == 1
                    st = 0;
                else
                    st = ends(idx-1);
                end
                % Select the elements of x that are effectively used in f
                xidx = (st+1):ends(idx);
                x = X(xidx,:);
                
                % X_a
                if j <= m
                    % 1=x_a, 2=y_a, 3=x_i {, y_i}
                    % mu(1)*xi.*ya - mu(3)*xa + rb;
                    fj = mu(1,:).*x(3,:).*x(2,:) - mu(3,:).*x(1,:);
                    % Boundary condition
                    if j == m
                        fj = fj + x(3,:).*mu(9,:)/s.hs;
                    end
                    
                % Y_a
                elseif m < j && j <= 2*m
                    % 1=x_a, 2=y_a {, x_i}, 3=y_i
                    % mu(2)*yi.*xan - mu(4)*ya;
                    fj = mu(2,:).*x(3,:).*x(1,:).^n - mu(4,:).*x(2,:);
                    
                % X_i
                elseif 2*m < j && j <= 3*m
                    % {x_a,} 1=y_a, 2=x_i {, y_i}
                    % -mu(1)*xi.*ya - mu(5)*xi + mu(7) - rb;
                    fj = -mu(1,:).*x(2,:).*x(1,:) - mu(5,:).*x(2,:) + mu(7,:);
                    % Boundary condition
                    if j == 3*m
                        fj = fj - x(2,:).*mu(9,:)/s.hs;
                    end
                    
                % Y_i
                else
                    % 1=x_a, {y_a, x_i}, 2=y_i
                    % -mu(2)*yi.*xan - mu(6)*yi + mu(8);
                    fj = -mu(2,:).*x(2,:).*x(1,:).^n - mu(6,:).*x(2,:) + mu(8,:);
                end
                fxj(idx,:) = fj;
            end
        end
    end
end

