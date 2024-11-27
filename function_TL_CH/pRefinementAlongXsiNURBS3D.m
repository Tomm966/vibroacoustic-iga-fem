function [Bnew,uKnotnew,weights_new,noPtsX_new,pb,T] = pRefinementAlongXsiNURBS3D(uKnot,B,p,ele_ord_xsi,noPtsX,noPtsY,noPtsZ,weights)
    % WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT
% uKnot: Knot vector in the ξ (u) direction.
% B: Matrix of control points in 3D space.
% p: Current degree in the ξ direction.
% ele_ord_xsi: Number of degrees to elevate in the ξ direction.
% noPtsX: Number of control points in the x direction.
% noPtsY: Number of control points in the y direction.
% noPtsZ: Number of control points in the z direction.
% weights: Weights associated with the control points.
% OUTPUT
% Bnew: Updated control points matrix after degree elevation.
% uKnotnew: New knot vector in the ξ direction after elevation.
% weights_new: New weights for the control points after degree elevation.
% noPtsX_new: Number of control points in the x direction after degree elevation.
% pb: New degree in the ξ direction after elevation.

    Bxsi = B([1:noPtsX],:);
    weightsxsi = weights([1:noPtsX],:);
    [~,uKnotnew,~,Txsi,pb] = pRefinement(uKnot,Bxsi,p,ele_ord_xsi,weightsxsi);
    noPtsX_new = length(uKnotnew)-pb-1;
    Bnew = zeros(noPtsX_new*noPtsY*noPtsZ,3);
    weights_new = zeros(noPtsX_new*noPtsY*noPtsZ,1);
    T    = zeros(noPtsX_new*noPtsY*noPtsZ,noPtsX*noPtsY*noPtsZ);
    for i = 1:noPtsY*noPtsZ
        numline = (i-1)*noPtsX+1: noPtsX*i;
        numline_new =  (i-1)*noPtsX_new+1:noPtsX_new*i;
        T(numline_new,numline)=Txsi;
    end
    weights_new = T*weights;
    Bnew = T*(B.*weights);
    Bnew = Bnew ./ weights_new;
end
