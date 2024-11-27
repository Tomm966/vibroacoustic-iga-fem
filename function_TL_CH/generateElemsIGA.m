function [nodesTot,elemsTot]=generateElemsIGA(NURBS)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% 
% The function takes a single input, NURBS, which is an object containing:
% controlPts: Control points defining the NURBS surface.
% weights: Weights associated with the control points.
% uknot, vknot, wknot: Knot vectors in each direction.
% p, q, r: Degrees of the NURBS in the respective directions.
% noPtsX, noPtsY, noPtsZ: Number of control points in each direction.
% OUTPUT:
% nodesTot: A matrix containing all physical nodes in 3D space.
% elemsTot: A matrix containing connectivity information for each element, defining which nodes make up the elements.

    
    noPtsX_1 = NURBS.noPtsX;
    noPtsY_1 = NURBS.noPtsY;
    noPtsZ_1 = NURBS.noPtsZ;
    controlPts = NURBS.controlPts;
    weights = NURBS.weights;
    uKnot_1 = NURBS.uknot;
    vKnot_1 = NURBS.vknot;
    wKnot_1 = NURBS.wknot;
    p = NURBS.p;
    q = NURBS.q;
    r = NURBS.r;


    uniqueUKnots   = unique(uKnot_1);
    uniqueVKnots   = unique(vKnot_1);
    uniqueWKnots   = unique(wKnot_1);
    noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
    noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.
    noElemsW       = length(uniqueWKnots)-1; % # of elements zeta dir.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % chan =
    %        1 2 3 4
    %        5 6 7 8
    %        9 10 11 12
    %        13 14 15 16
    % for a 4x2x2 control points
    
    
    chan  = zeros(noPtsZ_1,noPtsY_1,noPtsX_1);
    
    count = 1;
    
    for i=1:noPtsZ_1
        for j=1:noPtsY_1
            for k=1:noPtsX_1
                chan(i,j,k) = count;
                count       = count + 1;
            end
        end
    end
    
    
    % determine our element ranges and the corresponding
    % knot indices along each direction
    
    [elRangeU,elConnU] = buildConnectivity(p,uKnot_1,noElemsU);
    [elRangeV,elConnV] = buildConnectivity(q,vKnot_1,noElemsV);
    [elRangeW,elConnW] = buildConnectivity(r,wKnot_1,noElemsW);
    noElems = noElemsU * noElemsV * noElemsW;
    element = zeros(noElems,(p+1)*(q+1)*(r+1));
    e = 1;
    for w=1:noElemsW
        wConn = elConnW(w,:);
        for v=1:noElemsV
            vConn = elConnV(v,:);
            for u=1:noElemsU
                c = 1;
                uConn = elConnU(u,:);
                for i=1:length(wConn)
                    for j=1:length(vConn)
                        for k=1:length(uConn)
                            element(e,c) = chan(wConn(i),vConn(j),uConn(k));
                            c = c + 1;
                        end
                    end
                end
                e = e + 1;
            end
        end
    end
    index = zeros(noElems,3);
    count = 1;
    for i=1:size(elRangeW,1)
        for j=1:size(elRangeV,1)
            for k=1:size(elRangeU,1)
                index(count,1) = k;
                index(count,2) = j;
                index(count,3) = i;
                count = count + 1;
            end
        end
    end
    elemsTot = [];
    nodesTot = [];
    for e=1:noElems
        idu    = index(e,1);
        idv    = index(e,2);
        idw    = index(e,3);
        
        xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
        etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
        zetaE  = elRangeW(idw,:); % [zeta_k,zeta_k+1]
    
        nnl = 10;
        nodes1  = [linspace(xiE(1),xiE(2),nnl)' , ones(nnl,1)*etaE(1)            , ones(nnl,1)*zetaE(1)];
        nodes2  = [ones(nnl,1)*xiE(2)           , linspace(etaE(1),etaE(2),nnl)' , ones(nnl,1)*zetaE(1)];
        nodes3  = [linspace(xiE(1),xiE(2),nnl)' , ones(nnl,1)*etaE(2)            , ones(nnl,1)*zetaE(1)];
        nodes4  = [ones(nnl,1)*xiE(1)           , linspace(etaE(1),etaE(2),nnl)' , ones(nnl,1)*zetaE(1)];
        nodes5  = [ones(nnl,1)*xiE(1)           , ones(nnl,1)*etaE(1)            , linspace(zetaE(1),zetaE(2),nnl)'];
        nodes6  = [ones(nnl,1)*xiE(2)           , ones(nnl,1)*etaE(1)            , linspace(zetaE(1),zetaE(2),nnl)'];
        nodes7  = [ones(nnl,1)*xiE(2)           , ones(nnl,1)*etaE(2)            , linspace(zetaE(1),zetaE(2),nnl)'];
        nodes8  = [ones(nnl,1)*xiE(1)           , ones(nnl,1)*etaE(2)            , linspace(zetaE(1),zetaE(2),nnl)'];
        nodes9  = [linspace(xiE(1),xiE(2),nnl)' , ones(nnl,1)*etaE(1)            , ones(nnl,1)*zetaE(2)];
        nodes10 = [ones(nnl,1)*xiE(2)           , linspace(etaE(1),etaE(2),nnl)' , ones(nnl,1)*zetaE(2)];
        nodes11 = [linspace(xiE(1),xiE(2),nnl)' , ones(nnl,1)*etaE(2)            , ones(nnl,1)*zetaE(2)];
        nodes12 = [ones(nnl,1)*xiE(1)           , linspace(etaE(1),etaE(2),nnl)' , ones(nnl,1)*zetaE(2)];
        nodes = [nodes1;nodes2;nodes3;nodes4;nodes5;nodes6;nodes7;nodes8;nodes9;nodes10;nodes11;nodes12];
        elems = zeros(12*(nnl-1),2);
        for i = 1:12
            elemsi = [1:nnl-1 ; 2:nnl]';
            elems([1:nnl-1] + (i-1)*(nnl -1),:) = elemsi + (i-1)*nnl;
        end
        projcoord = nurb2proj(noPtsX_1*noPtsY_1*noPtsZ_1, controlPts, weights);
        dim = 4;
        nodes_physic = zeros(size(nodes,1),dim);
        for i = 1:size(nodes,1)
            nodes_physic(i,:) = SolidPoint(noPtsX_1-1,p, uKnot_1, ... 
                                           noPtsY_1-1,q, vKnot_1, ...
                                           noPtsZ_1-1,r, wKnot_1, ...
                                           projcoord,dim, ...,
                                           nodes(i,1),nodes(i,2),nodes(i,3));
        end
        nodes_physic(:,[1,2,3]) = nodes_physic(:,[1,2,3])./nodes_physic(:,4);
        elemsTot = [elemsTot ; elems + size(nodesTot,1)];
        nodesTot = [nodesTot ; nodes_physic];
        %matlab2VTK_net(nodes_physic(:,[1,2,3]),elems,zeros(3*size(nodes,1),1),ones(size(elems,1),1),['elem',num2str(e)])
    end
end