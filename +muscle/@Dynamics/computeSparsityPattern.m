function [SPK, SPalpha, SPLamDot] = computeSparsityPattern(this)
    % Computes all sorts of patterns simultaneously
    sys = this.fsys;
    mc = sys.Model.Config;
    fe_pos = mc.FEM;
    geo = fe_pos.Geometry;
    fe_press = mc.PressureFEM;
    pgeo = fe_press.Geometry;

    N = geo.NumNodes;
    
    % Jalpha
    if this.nfibres > 0
        iS = [];
        jS = [];
        columns_sarco_link = 53:56:56*this.nfibres;
    end
    
    % JLamDot
    if ~isempty(this.lambda_dot_pos)
        ildot = [];
        jldot = [];
    end
    
    num_elements = geo.NumElements;
    num_gausspoints = fe_pos.GaussPointsPerElem;
    dofs_pos = N*3;
    dofsperelem_displ = geo.DofsPerElement;
    dofsperelem_press = pgeo.DofsPerElement;
    
    % Compute indices vectors size for speed
    isize = num_elements*num_gausspoints*(...
            dofsperelem_displ*...
                (3*(3*dofsperelem_displ + dofsperelem_press)));
    i = zeros(isize,1,'int32');
    j = zeros(isize,1,'int32');    

    curoff = 0;
    for m = 1:num_elements
        elemidx_u = sys.idx_u_elems_local(:,:,m);
        elemidx_v = elemidx_u + dofs_pos;
        elemidx_p_glob = sys.idx_p_elems_local(:,m)+2*dofs_pos;
        inew = elemidx_u(:);
        one = ones(size(inew),'int32');
        for gp = 1:num_gausspoints
            for k = 1:dofsperelem_displ
                %% Grad_u K(u,v,w)
                step = 3*3*dofsperelem_displ;
                % x,y,zdim
                i(curoff+(1:step)) = [inew; inew; inew];
                j(curoff+(1:step)) = [one*elemidx_u(1,k)
                                      one*elemidx_u(2,k)
                                      one*elemidx_u(3,k)];
                curoff = curoff + step;
                
                %% Grad_v K(u,v,w) FIXME IF UNCOMMENTED
%                 if visc > 0
%                     % xdim
%                     i = [i; inew]; %#ok<*AGROW>
%                     j = [j; one*elemidx_velo(1,k)];
%                     % ydim
%                     i = [i; inew]; 
%                     j = [j; one*elemidx_velo(2,k)]; 
%                     % zdim
%                     i = [i; inew]; 
%                     j = [j; one*elemidx_velo(3,k)]; 
%                 end
            end
            
            %% Grad_w K(u,v,w)
            inew = elemidx_u(:);
            step = 3*dofsperelem_displ;
            for k = 1:dofsperelem_press
                i(curoff + (1:step)) = inew;
                j(curoff + (1:step)) = one*elemidx_p_glob(k); 
                curoff = curoff + step;
            end
            
            %% Jalpha pattern
            for k = 1:this.nfibres
                iS = [iS; elemidx_v(:)-dofs_pos];
                jS = [jS; ones(3*dofsperelem_displ,1)*columns_sarco_link(k)];
            end
            
            %% Check if change rate of lambda at a certain point should be tracked
            if ~isempty(this.lambda_dot_pos)
                k = find(this.lambda_dot_pos(1,:) == m & this.lambda_dot_pos(2,:) == gp);
                if ~isempty(k)
                    ildot = [ildot; k*ones(6*dofsperelem_displ,1)];
                    jldot = [jldot; elemidx_u(:); elemidx_v(:)];
                end
            end
        end
    end
    SPK = sparse(double(i),double(j),ones(size(i)),3*N,6*N+pgeo.NumNodes);
    
    % Remove values at dirichlet nodes
    SPK = SPK(1:3*N,:);
    SPK(:,sys.idx_uv_bc_glob) = [];
    SPK([sys.idx_u_bc_local; sys.idx_expl_v_bc_local],:) = [];
    SPK = logical(SPK);

    SPalpha = [];
    if this.nfibres > 0
        SPalpha = sparse(double(iS),double(jS),ones(size(iS)),3*N,this.nfibres*56);
        % Remove those that are connected to dirichlet values
        SPalpha([sys.idx_u_bc_local sys.idx_v_bc_local],:) = [];
        SPalpha = logical(SPalpha);
    end
    
    SPLamDot = [];
    if ~isempty(this.lambda_dot_pos)
        SPLamDot = sparse(ildot,double(jldot),true(size(ildot)),this.nfibres,6*N);
        SPLamDot(:,sys.idx_uv_bc_glob) = [];
    end
end

