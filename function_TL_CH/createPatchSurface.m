function [nodes_sur] = createPatchSurface(NURBS2D,nodes_param)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% NURBS2D: A structure containing the following fields related to the NURBS surface:
% noPtsX: The number of control points in the X direction.
% noPtsY: The number of control points in the Y direction.
% controlPts: An array of control points defining the NURBS surface, typically of size (noPtsX×noPtsY)×3 (or 4 if using homogeneous coordinates).
% weights: A vector of weights associated with the control points, with a length equal to the number of control points.
% uknot: The knot vector in the U direction, a sorted array of non-decreasing values.
% vknot: The knot vector in the V direction, similar to uknot.
% p: The degree of the NURBS surface in the U direction.
% q: The degree of the NURBS surface in the V direction.
% nodes_param: A matrix of size N×2, where each row contains the parametric coordinates (u,v) of a point on the NURBS surface.
% 
% OUTPUT:
% nodes_sur: A matrix of size N×4 containing the physical coordinates of the surface points corresponding to the parametric coordinates in nodes_param. Each row contains the X,Y,Z coordinates and a homogeneous coordinate.

    dim = 4;
    noPtsX = NURBS2D.noPtsX;
    noPtsY = NURBS2D.noPtsY;
    controlPts = NURBS2D.controlPts;
    weights = NURBS2D.weights;
    uKnot = NURBS2D.uknot;
    vKnot = NURBS2D.vknot;
    p = NURBS2D.p;
    q = NURBS2D.q;
    nodes_sur = zeros(size(nodes_param,1),dim);
    projcoord = nurb2proj(noPtsX*noPtsY, controlPts, weights);
    for j = 1:size(nodes_param,1)
       nodes_sur(j,:) = SurfacePoint(noPtsX-1,p, uKnot, ... 
                                     noPtsY-1,q, vKnot, ...
                                     projcoord,dim, ...
                                     nodes_param(j,1),nodes_param(j,2));
    end
    nodes_sur(:,[1,2,3]) = nodes_sur(:,[1,2,3])./nodes_sur(:,4); 
end