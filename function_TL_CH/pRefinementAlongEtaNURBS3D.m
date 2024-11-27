function [Bnew,wKnotnew,weights_new,noPtsZ_new,rb,T] = pRefinementAlongEtaNURBS3D(wKnot,B,r,ele_ord_eta,noPtsX,noPtsY,noPtsZ,weights)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% wKnot: Knot vector in the η (v) direction.
% B: Matrix of control points in 3D space.
% r: Current degree in the η direction.
% ele_ord_eta: Number of degrees to elevate in the η direction.
% noPtsX: Number of control points in the x direction.
% noPtsY: Number of control points in the y direction.
% noPtsZ: Number of control points in the z direction.
% weights: Weights associated with the control points.
% OUTPUT
% Bnew: Updated control points matrix after degree elevation.
% wKnotnew: New knot vector in the η direction after elevation.
% weights_new: New weights for the control points after degree elevation.
% noPtsZ_new: Number of control points in the z direction after degree elevation.
% rb: New degree in the η direction after elevation.

    Beta = B([1:noPtsX*noPtsY:((noPtsX*noPtsY*noPtsZ)-noPtsX*noPtsY+1)],:);
    weightseta = weights([1:noPtsX*noPtsY:((noPtsX*noPtsY*noPtsZ)-noPtsX*noPtsY+1)],:);
     [~,wKnotnew,~,Teta,rb] = pRefinement(wKnot,Beta,r,ele_ord_eta,weightseta);
    noPtsZ_new = length(wKnotnew)-rb-1;
    Bnew = zeros(noPtsX*noPtsZ_new*noPtsY,3);
    weights_new = zeros(noPtsX*noPtsY*noPtsZ_new,1);
    T    = zeros(noPtsX*noPtsZ_new*noPtsY,noPtsX*noPtsY*noPtsZ);
    for i = 1:noPtsX*noPtsY
        numline = (1:noPtsX*noPtsY:((noPtsX*noPtsY*noPtsZ)-noPtsX*noPtsY+1)) + (i-1);
        numline_new = (1:noPtsX*noPtsY:((noPtsX*noPtsZ_new*noPtsY)-noPtsX*noPtsY+1)) + (i-1);
        T(numline_new,numline)=Teta;
    end
    weights_new = T*weights;
    Bnew = T*(B.*weights);
    Bnew = Bnew ./ weights_new;
end

