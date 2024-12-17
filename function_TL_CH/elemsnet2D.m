function [elemsnet] = elemsnet2D(noPtsX,noPtsY)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% noPtsX: The number of control points in the X direction.
% noPtsY: The number of control points in the Y direction.
% 
% OUTPUT:
% net mesh


    elemsnet = zeros((noPtsX-1)*noPtsY + noPtsX*(noPtsY-1),2);
    if (noPtsX == 2 && noPtsY > 2)
        elemsnet([1:noPtsY],1) = [1:2:(2*noPtsY-1)]';
        elemsnet([1:noPtsY],2) = [2:2:(2*noPtsY)]';
        elemsnet([1:noPtsY-1]+noPtsY,1) = [1:2:2*noPtsY-3]';
        elemsnet([1:noPtsY-1]+noPtsY,2) = [3:2:2*noPtsY-1]';
        elemsnet([1:noPtsY-1]+2*noPtsY-1,1) = [1:2:2*noPtsY-3]'+1;
        elemsnet([1:noPtsY-1]+2*noPtsY-1,2) = [3:2:2*noPtsY-1]'+1;
    elseif (noPtsY == 2 && noPtsX > 2)
        %elemsnet = zeros(noPtsX+noPtsX-1,2);
        elemsnet([1:noPtsX],1) = [1:noPtsX]';
        elemsnet([1:noPtsX],2) = [1:noPtsX]'+noPtsX;
        elemsnet([1:noPtsX-1]+noPtsX,1) = [1:noPtsX-1]';
        elemsnet([1:noPtsX-1]+noPtsX,2) = [1:noPtsX-1]'+1;
        elemsnet([1:noPtsX-1]+2*noPtsX-1,1) = [1:noPtsX-1]'+noPtsX;
        elemsnet([1:noPtsX-1]+2*noPtsX-1,2) = [1:noPtsX-1]'+noPtsX+1;
    else
        for i = 1:noPtsY
            elemsnet( [1:noPtsX-1] + (i-1)*(noPtsY-1),1 ) = [1:noPtsX-1]' + (i-1)*noPtsY;   
            elemsnet( [1:noPtsX-1] + (i-1)*(noPtsY-1),2 ) = [1:noPtsX-1]' + (i-1)*noPtsY +1  ;
        end
        for i = 1:noPtsX  
            elemsnet([1:noPtsY-1] + (i-1)*(noPtsX-1) + (noPtsX-1)*noPtsY,1 ) = [[1:noPtsX:(noPtsX*noPtsY-noPtsX)] + i-1 ]'; 
            elemsnet([1:noPtsY-1] + (i-1)*(noPtsX-1) + (noPtsX-1)*noPtsY,2 ) = [[1:noPtsX:(noPtsX*noPtsY-noPtsX)] + i-1 + noPtsX ]'; 
        end
    end
end