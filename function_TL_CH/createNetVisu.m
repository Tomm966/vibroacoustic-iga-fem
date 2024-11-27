function [elemsnet] = createNet(noPtsX,noPtsY,noPtsZ)
 % WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% noPtsX: Number of points in the X-direction.
% noPtsY: Number of points in the Y-direction.
% noPtsZ: Number of points in the Z-direction.
% OUTPUT:
% elemsnet: A matrix representing the connectivity of the points in the X, Y, and Z directions, forming a network of elements.

    listnumX = [];
    for i = 1:noPtsY*noPtsZ
        listnumX = [ listnumX  ; ((i-1)*noPtsX+1) : i*noPtsX];
    end
    
    listnumZ = [];
    for i = 1:noPtsX*noPtsY
        listnumZcurrent = [];
        for j = 1:noPtsZ
            listnumZcurrent = [listnumZcurrent, i+(j-1)*noPtsX*noPtsY];
        end
        listnumZ = [listnumZ ; listnumZcurrent ];
    end
    
    listnumY = [];
    for k = 1:noPtsZ
        for i = 1:noPtsX %noPtsX*noPtsZ
            listnumYcurrent = [];
            for j = 1:noPtsY
                listnumYcurrent = [listnumYcurrent, i + (j-1)*noPtsX + (k-1)*noPtsX*noPtsY] ;
            end
            listnumY = [listnumY ; listnumYcurrent];
        end
    end

    elemsX = [];
    for i = 1:size(listnumX,1)
        for j = 1:size(listnumX,2)-1
            elemsX = [elemsX;listnumX(i,j),listnumX(i,j+1)];
        end
    end
    elemsY = [];
    for i = 1:size(listnumY,1)
        for j = 1:size(listnumY,2)-1
            elemsY = [elemsY;listnumY(i,j),listnumY(i,j+1)];
        end
    end
    elemsZ = [];
    for i = 1:size(listnumZ,1)
        for j = 1:size(listnumZ,2)-1
            elemsZ = [elemsZ;listnumZ(i,j),listnumZ(i,j+1)];
        end
    end
    elemsnet = [elemsX;elemsY;elemsZ];

end