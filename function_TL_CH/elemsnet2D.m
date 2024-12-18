function [elemsnet] = elemsnet2D(noPtsX,noPtsY)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% noPtsX: The number of control points in the X direction.
% noPtsY: The number of control points in the Y direction.
% 
% OUTPUT:
% net mesh connectivity

    elemsnet = [];
    %the dimension of elemsnet must be elemsnet= zeros((noPtsX-1)*noPtsY+noPtsX*(noPtsY-1),2)
    % Horizontal connectivity (along x)
    for row = 0:(noPtsY - 1)
        for col = 1:(noPtsX - 1)
            node1 = row * noPtsX + col;
            node2 = node1 + 1;
            elemsnet = [elemsnet; node1, node2];
        end
    end
    % Vertical connectivity (along y)
    for col = 1:noPtsX
        for row = 0:(noPtsY - 2)
            node1 = row * noPtsX + col;
            node2 = node1 + noPtsX;
            elemsnet = [elemsnet; node1, node2];
        end
    end
end