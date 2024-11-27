function [Bnew,vKnot_new,weights_new,noPtsY_new,T]=knotInsersionAlongNiNURBS3D(vKnot,ni_new,B,q,noPtsX,noPtsY,noPtsZ,weights)
    % WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT
% vKnot: The existing knot vector in the η direction.
% ni_new: The new knot value to be inserted.
% B: The current control points matrix (including weights).
% q: The degree of the basis functions in the η direction.
% noPtsX: The number of control points in the x direction.
% noPtsY: The number of control points in the y direction.
% noPtsZ: The number of control points in the z direction.
% weights: The weights associated with the control points.
% OUTPUT
% Bnew: The updated control points after knot insertion.
% vKnot_new: The updated knot vector after the insertion.
% weights_new: The updated weights for the new control points.
% noPtsY_new: The new number of control points in the y direction after insertion.
% T: Transformation matrix for mapping old control points to new control points.

noPtsY_new = noPtsY+1;
    B = B.*weights;
    Bni = B([1:noPtsX:noPtsX*noPtsY],:);
    % disp([1:3:noPtsY*noPtsX])
    % working in progress
    [~,Tni] = knotInserstion(vKnot,ni_new,Bni,q);
    Bnew = zeros(noPtsX*noPtsY_new*noPtsZ,3);
    weights_new = zeros(noPtsX*noPtsY_new*noPtsZ,1);
    T    = zeros(noPtsX*noPtsY_new*noPtsZ,noPtsX*noPtsY*noPtsZ);
    for j = 1:noPtsZ
        for i=1:noPtsX
        numline = (i:noPtsX:((noPtsX*noPtsY)-noPtsX+1+i))+(j-1)*noPtsX*noPtsY;
        numline_new = (1:noPtsX:((noPtsX*noPtsY_new)-noPtsX+1)) + (i-1)+(j-1)*noPtsX*noPtsY_new;
        %disp(numline);
        %disp(numline_new);
        weights_new(numline_new) =  Tni*weights(numline);
        Bnew(numline_new,:) = Tni*B(numline,:);
        T(numline_new,numline)=Tni;
        % keyboard
        end
     end
     Bnew= Bnew./weights_new;
     [vKnot_new] = newKnot(vKnot,ni_new);
end
