function Jg = getStateJacobianImpl(this, uvwdof, t)
%     J = this.getStateJacobianFD(uvwdof,t);
%     return;
    
    sys = this.fsys;
    
    mc = sys.Model.Config; 
    fe_pos = mc.FEM;
    geo = fe_pos.Geometry;
    fe_press = mc.PressureFEM;
    pgeo = fe_press.Geometry;

    N = geo.NumNodes;
    M = pgeo.NumNodes;
    fo = sys.NumFirstOrderDofs;
    dofsperelem_pos = geo.DofsPerElement;
    dofsperelem_press = pgeo.DofsPerElement;
    num_gausspoints = fe_pos.GaussPointsPerElem;
    num_elements = geo.NumElements;

    %% Precompute the size of i,j,s for speed
    relidx_press = 1:dofsperelem_press;
    totalsize = num_elements*num_gausspoints*dofsperelem_press;
    i = zeros(totalsize,1);
    j = zeros(totalsize,1);
    s = zeros(totalsize,1);
    
    idx_pos = sys.idx_u_elems_local;
    idx_press = sys.idx_p_elems_local;

    % Include dirichlet values to state vector
    uvwcomplete = sys.includeDirichletValues(t, uvwdof);
    
    cur_off = 0;
    for m = 1:num_elements
        elemidx_u_glob = idx_pos(:,:,m);
        elemidx_p = idx_press(:,m);
        elemidx_pressure3 = [elemidx_p
                             elemidx_p
                             elemidx_p];
       
        for gp = 1:num_gausspoints
            pos = 3*(gp-1)+1:3*gp;
            dtn = fe_pos.transgrad(:,pos,m);
            u = uvwcomplete(elemidx_u_glob);
            
            % Deformation gradient
            F = u * dtn;

            Finv = inv(F);

            %% Main node loop
            weight = fe_pos.GaussWeights(gp) * fe_pos.elem_detjac(m, gp);

            for k = 1:dofsperelem_pos
                % U_i^k = e_i dyad dPhik from script
                U1k = [dtn(k,:); 0 0 0; 0 0 0];
                U2k = [0 0 0; dtn(k,:); 0 0 0];
                U3k = [0 0 0; 0 0 0; dtn(k,:)];

                %% grad u g(u)
                precomp = weight * det(F) * fe_press.Ngp(:,gp,m);
                % Assign i index as whole for x,y,z (speed)
                i(cur_off + (1:3*dofsperelem_press)) = elemidx_pressure3;
                % dx
                j(cur_off + relidx_press) = elemidx_u_glob(1,k);
                s(cur_off + relidx_press) = sum(diag(Finv*U1k)) * precomp;%#ok
                cur_off = cur_off + dofsperelem_press;
                % dy
                j(cur_off + relidx_press) = elemidx_u_glob(2,k);
                s(cur_off + relidx_press) = sum(diag(Finv*U2k)) * precomp;%#ok
                cur_off = cur_off + dofsperelem_press;
                %dz
                j(cur_off + relidx_press) = elemidx_u_glob(3,k);
                s(cur_off + relidx_press) = sum(diag(Finv*U3k)) * precomp;%#ok
                cur_off = cur_off + dofsperelem_press;
            end
        end
    end
    Jg = sparse(i,j,s,M,6*N+M+fo);
    % Remove values at dirichlet nodes
    Jg(:,sys.idx_uv_bc_glob) = [];
end