function  [Bnew,uKnotnew,weightsnew,T,pb] = pRefinement(uKnot,B,p,n_elev,weights)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% uKnot: A vector containing the knot vector in the u-direction. This defines how the control points are associated with the parameter space.
% B: A matrix of control points. Each row represents a control point in the NURBS surface.
% p: The current polynomial degree of the NURBS in the u-direction.
% n_elev: The number of degrees to elevate.
% weights: A vector of weights associated with each control point.
% OUTPUT:
% Bnew: The updated matrix of control points after degree elevation.
% uKnotnew: The new knot vector after degree elevation.
% weightsnew: The updated weights after degree elevation.
% T: The transformation matrix used to compute the new control points.
% pb: The new degree after elevation.

    n  = length(uKnot)-p-1;
    pb = p+n_elev;
    uKnotcurrent = uKnot;
    for i = 1:pb-p
        uKnotnew = sort([unique(uKnotcurrent) , uKnotcurrent]);
        uKnotcurrent = uKnotnew;
    end
    uKnotnew;
    nb = length(uKnotnew)-pb-1;
    xsi_evaluation = linspace(0,1,nb);
    Mb = zeros(nb,nb);
    M  = zeros(nb,n);
    
    for j=1:nb-1
        knotSpanIndexnew = FindSpan(nb,pb,xsi_evaluation(j),uKnotnew);
        knotSpanIndex    = FindSpan(n ,p ,xsi_evaluation(j),uKnot);
        BasisFun(knotSpanIndexnew,xsi_evaluation(j),pb,uKnotnew);
        % one trick was here
        Mb(j,[knotSpanIndexnew-pb+1:knotSpanIndexnew+1]) = BasisFun(knotSpanIndexnew,xsi_evaluation(j),pb,uKnotnew)' ;
        M(j ,[knotSpanIndex-p+1:knotSpanIndex+1])  = BasisFun(knotSpanIndex,xsi_evaluation(j),p,uKnot)';
    end
    M(end,end) = 1;
    Mb(end,end) = 1;
%     tic
    T = sparse(inv(Mb)*M);
%     timeinverse = toc
    % the other trick is here
    Bnew = T*(B.*weights);
    weightsnew = T*weights;
    Bnew = Bnew ./ weightsnew;
end
