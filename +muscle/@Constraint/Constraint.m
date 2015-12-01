classdef Constraint < dscomponents.ACoreFun
    %CONSTRAINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ComputeUnassembled = false;
    end
    
    properties(Access=protected)
        % Reference to the full system
        fsys;
    end
    
    properties(SetAccess=private)
        idx_p_elems_unass;
        Sigma;
        fDim_unass;
    end
    
    methods
        function this = Constraint(sys)
            this = this@dscomponents.ACoreFun(sys);
        end
        
        function configUpdated(this)
            sys = this.fsys;
            mc = sys.Model.Config;
            if ~isempty(mc)
                this.fDim = sys.NumAlgebraicDofs;
                this.xDim = sys.NumAlgebraicDofs;
                this.JSparsityPattern = this.computeSparsityPattern;
                this.precomputeUnassembledData;
            end
        end
        
        function setSystem(this, sys)
            setSystem@dscomponents.ACoreFun(this, sys);
            if isa(sys,'models.ReducedSecondOrderSystem')
                this.fsys = sys.Model.FullModel.System;
            else
                this.fsys = sys;
            end
        end
        
        function evaluateCoreFun(~)
            error('Do not call. evaluate is implemented directly.');
        end
        
        function projected = project(this, V, W)
            projected = this.clone;
            projected = project@dscomponents.ACoreFun(this, V, W, projected);
        end
        
        function copy = clone(this)
            % Create new instance
            copy = models.muscle.Constraint(this.System);
            copy = clone@dscomponents.ACoreFun(this, copy);
            copy.ComputeUnassembled = this.ComputeUnassembled;
            copy.fsys = this.fsys;
            copy.idx_p_elems_unass = this.idx_p_elems_unass;
            copy.Sigma = this.Sigma;
            copy.fDim_unass = this.fDim_unass;
        end
    end
    
    methods(Access=private)
        function precomputeUnassembledData(this)
            sys = this.System;
            mc = sys.Model.Config;
            fem = mc.PressureFEM;
            geo = fem.Geometry;
            
            [I, ~] = find(fem.Sigma);
            
            n = numel(I);
            S = sparse(I,1:n,ones(n,1),geo.NumNodes,n);
            
            % Take out nodes with dirichlet BC on output side
%             S(sys.idx_v_bc_local,:) = [];
            % Find corresponding unassembled dofs that would be ignored
            % (due to dirichlet velocity values, pressure dirichlet not
            % implemented)
%             bc_unass = find(sum(S,1) == 0);
            % Remove them, too. The unassembled evaluation also removes the
            % corresponding entries of the unassembled vector.
%             S(:,bc_unass) = [];
            
            this.Sigma = S;
            
            this.fDim_unass = size(S,2);
            
            % Create boolean array that determines which unassembled dofs
            % belong to which element
            % dg dofs
            hlp = repmat(1:geo.NumElements,geo.DofsPerElement,1);
            hlp = hlp(:);
            ass = false(geo.NumElements,length(hlp));
            for k = 1:geo.NumElements
                ass(k,:) = hlp == k;
            end
            this.idx_p_elems_unass = ass;
        end
    end
    
end

