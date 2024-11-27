function [Bnew,uKnot_new,weights_new,noPtsX_new,T] = knotInsersionAlongXsiNURBS3D(uKnot,xsi_new,B,p,noPtsX,noPtsY,noPtsZ,weights)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU
    
% INPUT
% uKnot: The existing knot vector in the ξ direction.
% xsi_new: The new knot value to be inserted.
% B: The current control points matrix (including weights).
% p: The degree of the basis functions in the ξ direction.
% noPtsX: The number of control points in the x direction.
% noPtsY: The number of control points in the y direction.
% noPtsZ: The number of control points in the z direction.
% weights: The weights associated with the control points.
% OUTPUT
% Bnew: The updated control points after knot insertion.
% uKnot_new: The updated knot vector after the insertion.
% weights_new: The updated weights for the new control points.
% noPtsX_new: The new number of control points in the x direction after insertion.
% T: Transformation matrix for mapping old control points to new control points.

noPtsX_new = noPtsX+1;
    B = B.*weights;
    Bxsi = B([1:noPtsX],:);
    % disp([1:noPtsX])
    % working in progress
    [~,Txsi] = knotInserstion(uKnot,xsi_new,Bxsi,p);
    Bnew = zeros(noPtsX_new*noPtsY*noPtsZ,3);
    weights_new = zeros(noPtsX_new*noPtsY*noPtsZ,1);
    T    = zeros(noPtsX_new*noPtsY*noPtsZ,noPtsX*noPtsY*noPtsZ);
    for i = 1:noPtsY*noPtsZ
        numline = (i-1)*noPtsX+1: noPtsX*i;
        numline_new =  (i-1)*noPtsX_new+1:noPtsX_new*i;
        % disp(numline)
        % disp(numline_new)
        weights_new(numline_new) =  Txsi*weights(numline);
        Bnew(numline_new,:) = Txsi*B(numline,:);
        T(numline_new,numline)=Txsi;
    end
    Bnew= Bnew./weights_new;
    [uKnot_new] = newKnot(uKnot,xsi_new);
end