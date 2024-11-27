function [Bnew,vKnotnew,weights_new,noPtsY_new,qb,T] = pRefinementAlongNiNURBS3D(vKnot,B,q,ele_ord_ni,noPtsX,noPtsY,noPtsZ,weights)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT
% vKnot: Knot vector in the η (v) direction.
% B: Matrix of control points in 3D space.
% q: Current degree in the η direction.
% ele_ord_ni: Number of degrees to elevate in the η direction.
% noPtsX: Number of control points in the x direction.
% noPtsY: Number of control points in the y direction.
% noPtsZ: Number of control points in the z direction.
% weights: Weights associated with the control points.
% OUTPUT
% Bnew: Updated control points matrix after degree elevation.
% vKnotnew: New knot vector in the η direction after elevation.
% weights_new: New weights for the control points after degree elevation.
% noPtsY_new: Number of control points in the y direction after degree elevation.
% qb: New degree in the η direction after elevation.

    Bni = B([1:noPtsX:noPtsX*noPtsY],:);
    weightsni = weights([1:noPtsX:noPtsX*noPtsY],:);
    [~,vKnotnew,~,Tni,qb] = pRefinement(vKnot,Bni,q,ele_ord_ni,weightsni);
    noPtsY_new = length(vKnotnew)-qb-1;
    Bnew = zeros(noPtsX*noPtsY_new*noPtsZ,3);
    weights_new = zeros(noPtsX*noPtsY_new*noPtsZ,1);
    T    = zeros(noPtsX*noPtsY_new*noPtsZ,noPtsX*noPtsY*noPtsZ);
    for j = 1:noPtsZ
        for i=1:noPtsX
        numline = (i:noPtsX:((noPtsX*noPtsY)-noPtsX+1+i))+(j-1)*noPtsX*noPtsY;
        numline_new = (1:noPtsX:((noPtsX*noPtsY_new)-noPtsX+1)) + (i-1)+(j-1)*noPtsX*noPtsY_new;
        T(numline_new,numline)=Tni;
        end
    end
    weights_new = T*weights;
    Bnew = T*(B.*weights);
    Bnew = Bnew ./ weights_new;
end

