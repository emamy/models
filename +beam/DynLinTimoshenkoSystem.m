classdef DynLinTimoshenkoSystem < models.BaseDynSystem
% DynLinTimoshenkoSystem: 
%
%
%
% @author Daniel Wirtz @date 2011-09-20
%
% @new{0,5,dw,2011-09-20} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=private)
        % Mapping of effectively used points in output plotting (yet to be
        % moved elsewhere)
        KnotenIndex;
    end
    
    methods
        function this = DynLinTimoshenkoSystem(model)
            this = this@models.BaseDynSystem(model);
            dim = model.dim;
            this.f = models.beam.DynLinTimoshenkoCoreFun(model);
            
            this.MaxTimestep = [];
            
            %% Initial values
            x0 = 0*ones(dim,1);
            this.x0 = dscomponents.ConstInitialValue(x0);
            
            %% Input
            B = ones(dim,3);
            this.B = dscomponents.LinearInputConv(B);
            this.Inputs{1} = @(t)[sin(t); 0; 0];
            this.Inputs{2} = @(t)[sin(t); cos(t); 0];
            this.Inputs{3} = @(t)[0; t; cos(t)];
            
            %% Output setup (dim = 3*space + 3*velo + 1*temp)
            %xdim = (dim - (dim/7))/2;
            d = sparse(1:dim,1:dim,1);
            d(repmat(logical([0 0 0 1 1 1 0]'),dim/7,1),:) = [];
            %d = reshape(repmat([1 1 1 0 0 0 1]', dim/7,1),[],1);
            %d = [d; 0*d];
            this.C = dscomponents.LinearOutputConv(d);
        end
        
        function prepareConstants(this, mu, inputidx)
            % Hier KnotenIndex (== knoten_index) zuweisen (oder muss in TimoFun passieren?)
            this.KnotenIndex = 1:(this.Model.dim*4/7);
        end
    end
    
end