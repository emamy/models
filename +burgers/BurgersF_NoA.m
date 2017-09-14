classdef BurgersF_NoA < dscomponents.ACompEvalCoreFun
% BurgersF: 
%
%
%
% @author Daniel Wirtz @date 2012-04-24
%
% @new{0,6,dw,2012-04-24} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        Ax;
    end
    
    methods
        function this = BurgersF_NoA(sys)
            this = this@dscomponents.ACompEvalCoreFun(sys);
            this.CustomProjection = false;
            this.TimeDependent = false;
        end
        
        function fx = evaluateCoreFun(this, x, ~, ~)
            fx = - x.*(this.Ax*x);
        end
                
        function J = getStateJacobian(this, x, ~, ~)
            hlp = bsxfun(@times,this.Ax,x);
            J = -hlp - spdiags(this.Ax*x,0,size(x,1),size(x,1));
        end
        
        function newDim(this)
            m = this.System.Model;
            n = m.Dimension;
            dx = (m.Omega(2) - m.Omega(1))/(n+1);
            e = ones(n,1);
            d1 = e/(2*dx);
            this.Ax = spdiags([-d1 0*d1  d1], -1:1, n, n);
            this.JSparsityPattern = spdiags([e e  e], -1:1, n, n);
            spy(this.JSparsityPattern);
            this.xDim = n;
            this.fDim = n;
        end
        
%         function target = project(this, V, W)
%             target = this.clone;
%             target = project@dscomponents.ACoreFun()
%         end
        
        function copy = clone(this)
            copy = clone@dscomponents.ACompEvalCoreFun(this, models.burgers.BurgersF_NoA(this.System));
            copy.Ax = this.Ax;
        end
    end
    
    methods(Access=protected)
        function fxj = evaluateComponents(this, pts, ends, argidx, self, X, ~, ~)
            % Evaluates the burgers nonlinearity pointwise.
            %
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
            % Anzahl Punkte, die ausgewertet werden; Zeit
            fxj = zeros(length(pts),size(X,2));
            % Iteration über Anzahl der auszuwertenden f_i
            for idx=1:length(pts)
                %derzeitiger globaler Index von f_i, das ausgewertet wird
                pt = pts(idx);
                if idx == 1
                    st = 0;
                else
                    st = ends(idx-1);
                end
                % Select the elements of x that are effectively used in f
                % ends(i-1) + 1: ends(i+1)
                xidx = (st+1):ends(idx);
                % wirklich benötigte Werte von X für f_i über alle Zeitschritte
                x = X(xidx,:);
                %disp(self(xidx));
                % gebe nur punkt zurück, der auch teil des outputs von f
                % xidx ist auch nur lokaler Index
                %disp(x(self(xidx),:));
                % globaler Index
                %disp(argidx(xidx));
                % Teil der Matrix, der verwendet wird
                %disp(this.Ax(pt,argidx(xidx)));
                % matrix mit vielen punkten: wähle nur den teil, der die
                % involvierten punkte betrifft
                % disp(this.Ax(pt,argidx(xidx)));
                % fxj(idx,:) = nicht globale Position
                % einzige Möglichkeit
                fxj(idx,:) = - x(self(xidx),:) .* (this.Ax(pt,argidx(xidx))*x);
            end
            
        end
        
        function fxj = evaluateComponentsMulti(this, varargin)
            fxj = this.evaluateComponents(varargin{:});
        end
    end
end