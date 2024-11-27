function [elemsnet] = elemsnet2D(noPtsX,noPtsY)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% noPtsX: The number of control points in the X direction.
% noPtsY: The number of control points in the Y direction.
% OUTPUT:
% elemsnet: A matrix of size (noPtsX−1)×noPtsY+noPtsX×(noPtsY−1)×2 where each row represents an edge element in the network. Each element is defined by a pair of nodes (indices) forming the edges of a grid.

    elemsnet = zeros((noPtsX-1)*noPtsY + noPtsX*(noPtsY-1),2);
    for i = 1:noPtsY
        elemsnet( [1:noPtsX-1] + (i-1)*(noPtsY-1),1 ) = [1:noPtsX-1]' + (i-1)*noPtsY    ;
        elemsnet( [1:noPtsX-1] + (i-1)*(noPtsY-1),2 ) = [1:noPtsX-1]' + (i-1)*noPtsY +1  ;
    end
    
    for i = 1:noPtsX  
        elemsnet([1:noPtsY-1] + (i-1)*(noPtsX-1) + (noPtsX-1)*noPtsY,1 ) = [[1:noPtsX:(noPtsX*noPtsY-noPtsX)] + i-1 ]' ;
        elemsnet([1:noPtsY-1] + (i-1)*(noPtsX-1) + (noPtsX-1)*noPtsY,2 ) = [[1:noPtsX:(noPtsX*noPtsY-noPtsX)] + i-1 + noPtsX ]' ;
    end
end
