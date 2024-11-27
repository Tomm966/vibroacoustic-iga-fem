function [nodes,nodes_physic,elems] = createPatchEdges(NURBS,nnl)
    % WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% NURBS: A structure containing NURBS-related information with the following fields:
% noPtsX: Number of points in the X-direction.
% noPtsY: Number of points in the Y-direction.
% noPtsZ: Number of points in the Z-direction.
% controlPts: Control points of the NURBS surface.
% weights: Weights associated with the control points.
% uknot: Knot vector in the U-direction.
% vknot: Knot vector in the V-direction.
% wknot: Knot vector in the W-direction.
% p: Degree in the U-direction.
% q: Degree in the V-direction.
% r: Degree in the W-direction.
% nnl: Number of nodes along the edges of the NURBS surface.

% OUTPUT:
% nodes: A matrix containing the coordinates of the nodes defined on the edges of the NURBS surface.
% nodes_physic: A matrix containing the physical coordinates of the nodes after projecting them through the NURBS representation.
% elems: A matrix representing the connectivity of the nodes, indicating which nodes are connected to form line elements along the edges of the NURBS surface.

    %----------------------------------
    % Extract data from NURBS structure
    noPtsX = NURBS.noPtsX;
    noPtsY = NURBS.noPtsY;
    noPtsZ = NURBS.noPtsZ;
    controlPts = NURBS.controlPts;
    weights = NURBS.weights;
    uKnot = NURBS.uknot;
    vKnot = NURBS.vknot;
    wKnot = NURBS.wknot;
    p = NURBS.p;
    q = NURBS.q;
    r = NURBS.r;


    nodes1 = [linspace(0,1,nnl)' , zeros(nnl,1)       , zeros(nnl,1)];
    nodes2 = [ones(nnl,1)        , linspace(0,1,nnl)' , zeros(nnl,1)];
    nodes3 = [linspace(0,1,nnl)' , ones(nnl,1)        , zeros(nnl,1)];
    nodes4 = [zeros(nnl,1)       , linspace(0,1,nnl)' , zeros(nnl,1)];
    nodes5 = [zeros(nnl,1)       , zeros(nnl,1)       , linspace(0,1,nnl)'];
    nodes6 = [ones(nnl,1)        , zeros(nnl,1)       , linspace(0,1,nnl)'];
    nodes7 = [ones(nnl,1)        , ones(nnl,1)        , linspace(0,1,nnl)'];
    nodes8 = [zeros(nnl,1)       , ones(nnl,1)        , linspace(0,1,nnl)'];
    nodes9  = [linspace(0,1,nnl)' , zeros(nnl,1)       , ones(nnl,1)];
    nodes10 = [ones(nnl,1)        , linspace(0,1,nnl)' , ones(nnl,1)];
    nodes11 = [linspace(0,1,nnl)' , ones(nnl,1)        , ones(nnl,1)];
    nodes12 = [zeros(nnl,1)       , linspace(0,1,nnl)' , ones(nnl,1)];
    nodes = [nodes1;nodes2;nodes3;nodes4;nodes5;nodes6;nodes7;nodes8;nodes9;nodes10;nodes11;nodes12];
    elems = zeros(12*(nnl-1),2);
    for i = 1:12
        elemsi = [1:nnl-1 ; 2:nnl]';
        elems([1:nnl-1] + (i-1)*(nnl -1),:) = elemsi + (i-1)*nnl;
    end
    projcoord = nurb2proj(noPtsX*noPtsY*noPtsZ, controlPts, weights);
    dim = 4;
    nodes_physic = zeros(size(nodes,1),dim);
    for i = 1:size(nodes,1)
        nodes_physic(i,:) = SolidPoint(noPtsX-1,p, uKnot, ... 
                                       noPtsY-1,q, vKnot, ...
                                       noPtsZ-1,r, wKnot, ...
                                       projcoord,dim, ...,
                                       nodes(i,1),nodes(i,2),nodes(i,3));
    end
    nodes_physic(:,[1,2,3]) = nodes_physic(:,[1,2,3])./nodes_physic(:,4);

end