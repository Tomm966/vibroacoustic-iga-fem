function [field_vol] = createPatchField(NURBS,nodes_param,Ux,Uy,Uz)
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
% nodes_param: A matrix of size N×3 containing the parametric coordinates of N nodes on the NURBS surface.
% Ux: A vector representing the X-component of a field defined at the control points.
% Uy: A vector representing the Y-component of a field defined at the control points.
% Uz: A vector representing the Z-component of a field defined at the control points.
% OUTPUT:
% field_vol: A matrix of size N×4 where each row corresponds to the physical coordinates of the field defined at the nodes in nodes param. The four columns contain the X,Y,Z coordinates, and a homogeneous coordinate.

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
    field_vol = zeros(size(nodes_param,1),dim);
    projfield = nurb2proj(noPtsX*noPtsY*noPtsZ, [Ux,Uy,Uz], weights);
    for j = 1:size(nodes_param,1)
       field_vol(j,:) = SolidPoint(noPtsX-1,p, uKnot, ... 
                                   noPtsY-1,q, vKnot, ...
                                   noPtsZ-1,r, wKnot, ...
                                   projfield,dim, ...
                                   nodes_param(j,1),nodes_param(j,2),nodes_param(j,3));
    end
    field_vol(:,[1,2,3]) = field_vol(:,[1,2,3])./field_vol(:,4); 

end