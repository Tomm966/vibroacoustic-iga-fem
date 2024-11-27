function [B_new,T] = knotInserstion(uKnot,xsi_new,B,p)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT
% uKnot: The existing knot vector along the ξ direction.
% xsi_new: The new knot value to be inserted into the knot vector.
% B: The current control points matrix (including weights).
% p: The degree of the basis functions in the ξ direction.
% OUTPUT
% B_new: The updated control points after knot insertion.
% T: The transformation matrix that maps the old control points to the new control points.

    n = length(uKnot) - p - 1;
    indices = FindSpan(n ,p ,xsi_new,uKnot);
    k = indices-p+1;
    T = zeros(size(B,1)+1,size(B,1));
    for i=1:k
        T(i,i) = 1;
    end
    
    for i = k+p+1:(size(T,1))
        T(i,i-1) = 1;
    end
    
    for i= k+1:k+p
        ai = (xsi_new - uKnot(i) )/(uKnot(i+p) - uKnot(i));
        T(i,i-1)  = (1 - ai);
        T(i,i)    = ai;
    end
    B_new = T*B;
end