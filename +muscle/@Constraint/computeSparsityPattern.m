function SP = computeSparsityPattern(this)
    % Computes all sorts of patterns simultaneously
    sys = this.fsys;
    mc = sys.Model.Config;
    fe_pos = mc.FEM;
    geo = fe_pos.Geometry;
    fe_press = mc.PressureFEM;
    pgeo = fe_press.Geometry;

    num_elements = geo.NumElements;
    num_gausspoints = fe_pos.GaussPointsPerElem;
    dofsperelem_displ = geo.DofsPerElement;
    dofsperelem_press = pgeo.DofsPerElement;
    fo = sys.NumFirstOrderDofs;
    
    % Compute indices vectors size for speed
    itotal = num_elements*num_gausspoints...
        *(dofsperelem_press*(3*dofsperelem_displ)) ;
    i = zeros(itotal,1,'int32');
    j = zeros(itotal,1,'int32');    

    curoff = 0;
    pones = ones(dofsperelem_press,1,'int32');
    for m = 1:num_elements
        elemidx_u = sys.idx_u_elems_local(:,:,m);
        elemidx_p = sys.idx_p_elems_local(:,m);
        for gp = 1:num_gausspoints
            for k = 1:dofsperelem_displ
                %% grad u g(u)
                step = 3*dofsperelem_press;
                % dx,dy,dz
                i(curoff + (1:step)) = [elemidx_p(:)
                                        elemidx_p(:)
                                        elemidx_p(:)];
                j(curoff + (1:step)) = [pones*elemidx_u(1,k)
                                        pones*elemidx_u(2,k)
                                        pones*elemidx_u(3,k)];
                curoff = curoff + step;
            end
        end
    end
    % Create g constraint pattern
    M = pgeo.NumNodes;
    SP = sparse(double(i),double(j),ones(size(i)),M,6*geo.NumNodes+fo+M);
    % Reomve BC affected nodes
    SP(:,sys.idx_uv_bc_glob) = [];
    SP = logical(SP);
end

