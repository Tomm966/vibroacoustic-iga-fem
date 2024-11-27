function [nodes_vol] = createPatchVolume(NURBS,nodes_param)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% NURBS: A structure containing the following fields related to the NURBS volume:
% noPtsX: The number of control points in the X direction.
% noPtsY: The number of control points in the Y direction.
% noPtsZ: The number of control points in the Z direction.
% controlPts: An array of control points defining the NURBS volume, typically of size (noPtsX×noPtsY×noPtsZ)×3 (or 4 if using homogeneous coordinates).
% weights: A vector of weights associated with the control points, with a length equal to the number of control points.
% uknot: The knot vector in the U direction, a sorted array of non-decreasing values.
% vknot: The knot vector in the V direction, similar to uknot.
% wknot: The knot vector in the W direction, similar to uknot.
% p: The degree of the NURBS volume in the U direction.
% q: The degree of the NURBS volume in the V direction.
% r: The degree of the NURBS volume in the W direction.
% nodes_param: A matrix of size N×3, where each row contains the parametric coordinates (u,v,w) of a point in the NURBS volume.
% 
% OUTPUT:
% nodes_vol: A matrix of size N×4 containing the physical coordinates of the volume points corresponding to the parametric coordinates in nodes_param. Each row contains the X,Y,Z coordinates and a homogeneous coordinate.

    dim = 4;
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
    nodes_vol = zeros(size(nodes_param,1),dim);
    projcoord = nurb2proj(noPtsX*noPtsY*noPtsZ, controlPts, weights);
    for j = 1:size(nodes_param,1)
       nodes_vol(j,:) = SolidPoint(noPtsX-1,p, uKnot, ... 
                                   noPtsY-1,q, vKnot, ...
                                   noPtsZ-1,r, wKnot, ...
                                   projcoord,dim, ...
                                   nodes_param(j,1),nodes_param(j,2),nodes_param(j,3));
    end
    nodes_vol(:,[1,2,3]) = nodes_vol(:,[1,2,3])./nodes_vol(:,4); 

end