function [Bnew,wKnot_new,weights_new,noPtsZ_new,T]=knotInsersionAlongEtaNURBS3D(wKnot,eta_new,B,r,noPtsX,noPtsY,noPtsZ,weights)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT
% wKnot: The existing knot vector in the ζ direction.
% eta_new: The new knot value to be inserted.
% B: The current NURBS control points matrix (including weights).
% r: The degree of the basis functions in the η direction.
% noPtsX: The number of control points in the x direction.
% noPtsY: The number of control points in the y direction.
% noPtsZ: The number of control points in the z direction.
% weights: The weights associated with the control points.
% OUTPUT
% Bnew: The updated control points after knot insertion.
% wKnot_new: The updated knot vector after the insertion.
% weights_new: The updated weights for the new control points.
% noPtsZ_new: The new number of control points in the z direction after insertion.
% T: Transformation matrix for mapping old control points to new control points.   

noPtsZ_new = noPtsZ+1;
    B = B.*weights;
    Beta = B([1:noPtsX*noPtsY:((noPtsX*noPtsY*noPtsZ)-noPtsX*noPtsY+1)],:);
    % disp([1:noPtsX*noPtsY:((noPtsX*noPtsY*noPtsZ)-noPtsX*noPtsY+1)])
    % working in progress
    [~,Teta] = knotInserstion(wKnot,eta_new,Beta,r);
    Bnew = zeros(noPtsX*noPtsZ_new*noPtsY,3);
    weights_new = zeros(noPtsX*noPtsY*noPtsZ_new,1);
    T    = zeros(noPtsX*noPtsZ_new*noPtsY,noPtsX*noPtsY*noPtsZ);
    for i = 1:noPtsX*noPtsY
        %numline = (1:noPtsX:(noPtsX-1)*noPtsY+1) + (i-1)
        %numline_new = (1:noPtsX:(noPtsX-1)*(noPtsY_new+1)) + (i-1)
        numline = (1:noPtsX*noPtsY:((noPtsX*noPtsY*noPtsZ)-noPtsX*noPtsY+1)) + (i-1);
        numline_new = (1:noPtsX*noPtsY:((noPtsX*noPtsZ_new*noPtsY)-noPtsX*noPtsY+1)) + (i-1);
        % disp(numline)
        % disp(numline_new)
        weights_new(numline_new) =  Teta*weights(numline);
        Bnew(numline_new,:) = Teta*B(numline,:);
        T(numline_new,numline)=Teta;
    end
    Bnew= Bnew./weights_new;
    [wKnot_new] = newKnot(wKnot,eta_new);
end

