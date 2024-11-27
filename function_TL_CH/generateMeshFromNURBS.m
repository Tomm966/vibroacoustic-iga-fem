function [nodes_unique,elems_unique,elems_skin_unique] = generateMeshFromNURBS(nex,ney,nez,eletype,NURBS)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU
% INPUT:
% nex, ney, nez: Number of elements along the x, y, and z axes, respectively.
% eletype: Determines whether to use a linear mesh (eletype == 1) or a quadratic mesh (eletype != 1).
% NURBS: A structure that contains the NURBS geometry information, including control points, weights, and knot vectors.

% OUTPUT:
% nodes_unique: A list of unique node coordinates in the physical space.
% elems_unique: Element connectivity in terms of the unique nodes.
% elems_skin_unique: Surface element connectivity in terms of the unique nodes.

%====================================================================
    % Mesh the parametrix space with hexa20
    %====================================================================
    Lx = 1.0;
    Ly = 1.0;
    Lz = 1.0;
    %-------------------------------
    if eletype == 1
        [nodes_param,elems_param,elems_skin_param] = generateFEMeshLinear(Lx,Ly,Lz,nex,ney,nez);
    else
        [nodes_param,elems_param,elems_skin_param] = generateFEMeshQuadratique20(Lx,Ly,Lz,nex,ney,nez);
    end
    %
    
    %====================================================================
    % Mesh the physical space
    %====================================================================
    elems = [];
    nodes = [];
    elems_skin = [];
    for j = 1:length(NURBS)
        p = NURBS.p;
        q = NURBS.q;
        r = NURBS.r;
        noPtsX = NURBS(j).noPtsX;
        noPtsY = NURBS(j).noPtsY;
        noPtsZ = NURBS(j).noPtsZ;
        controlPts = NURBS(j).controlPts;
        weights = NURBS(j).weights;
        uKnot = NURBS(j).uknot;
        vKnot = NURBS(j).vknot;
        wKnot = NURBS(j).wknot;
        projcoord = nurb2proj(noPtsX*noPtsY*noPtsZ, controlPts, weights);
        dim = 4;
        nodes_physic = zeros(size(nodes_param,1),dim);
        for i = 1:size(nodes_param,1)
            nodes_physic(i,:) = SolidPoint(noPtsX-1,p, uKnot, ... 
                                           noPtsY-1,q, vKnot, ...
                                           noPtsZ-1,r, wKnot, ...
                                           projcoord,dim, ...,
                                           nodes_param(i,1),nodes_param(i,2),nodes_param(i,3));
        end
        nodes_physic(:,[1,2,3]) = nodes_physic(:,[1,2,3])./nodes_physic(:,4); 
        elems = [elems ; elems_param + size(nodes,1)];
        elems_skin = [elems_skin; elems_skin_param + size(nodes,1) ];
        nodes = [nodes ; nodes_physic];
    end
    [nodes_unique,~,renum] = uniquetol(nodes(:,[1,2,3]),'ByRows',1e-8);
    elems_unique = renum(elems);
    elems_skin_unique = renum(elems_skin);
end

