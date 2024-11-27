function [uKnot_new] = newKnot(uKnot,xsi_new)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT
% uKnot: A vector containing the existing knot values in non-decreasing order.
% xsi_new: A single scalar value representing the new knot to be inserted into the uKnot vector.
% OUTPUT
% uKnot_new: A new knot vector that includes the inserted knot while maintaining the non-decreasing order.
 
    indices = find(uKnot <= xsi_new);
    interval_start = indices(end);
    interval_end  = interval_start +1;
    uKnot_new = [uKnot(1:interval_start), ...
                 xsi_new, ...
                 uKnot(interval_end:end)];
end