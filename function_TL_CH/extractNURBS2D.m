function [NURBSwith2D] = extractNURBS2D(NURBS)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% NURBS: A structure array containing one or more 3D NURBS patches. Each entry should include the following fields:
% controlPts: Control points of the 3D NURBS.
% noPtsX: Number of control points in the X direction.
% noPtsY: Number of control points in the Y direction.
% noPtsZ: Number of control points in the Z direction.
% p: Degree in the X direction.
% q: Degree in the Y direction.
% r: Degree in the Z direction.
% uknot, vknot, wknot: Knot vectors in the X,Y, and Z directions respectively.
% weights: Weights associated with the control points.
% OUTPUT:
% NURBSwith2D: A modified structure array that includes an additional field, NURBS2D, for each 3D NURBS patch. This field contains an array of 2D NURBS representations for each of the six boundary faces of the 3D NURBS.

    for i = 1:length(NURBS)
        % Read informations for each patchs
        controlPts = NURBS(i).controlPts;
        noPtsX     = NURBS(i).noPtsX;
        noPtsY     = NURBS(i).noPtsY;
        noPtsZ     = NURBS(i).noPtsZ;
        p          = NURBS(i).p;
        q          = NURBS(i).q;
        r          = NURBS(i).r;
        uknot      = NURBS(i).uknot;
        vknot      = NURBS(i).vknot;
        wknot      = NURBS(i).wknot;
        weights    = NURBS(i).weights;
        
        % Creation of 6 2D NURBS for each 3D patchs
        NURBS2D    = [];
    
        % Face 1 : bottom zeta = 0 (external normal 1)
        NURBS2D(1).noPtsX = noPtsX;
        NURBS2D(1).noPtsY = noPtsY;
        NURBS2D(1).p      = p;
        NURBS2D(1).q      = q;
        B = zeros(noPtsX*noPtsY,3);
        for k = 1:noPtsY
           B([1:noPtsX]' + (k-1)*noPtsX , :) = controlPts([[1+noPtsX-1:-1:1] + (k-1)*noPtsX]',:);
        end
        NURBS2D(1).controlPts = B;
        NURBS2D(1).elemsnet = elemsnet2D(NURBS2D(1).noPtsX,NURBS2D(1).noPtsY);
        NURBS2D(1).PID    = 1;
        NURBS2D(1).uknot  = sort(1 - uknot); 
        NURBS2D(1).vknot  = vknot;
        w = zeros(noPtsX*noPtsY,1);
        for k = 1:noPtsY
            w([1:noPtsX]' + (k-1)*noPtsX , :) = weights([[1+noPtsX-1:-1:1] + (k-1)*noPtsX]',:);
        end
        NURBS2D(1).weights= w;
    
    
        % Face 2 : up zeta = 1 (external normal 1)
        NURBS2D(2).noPtsX = noPtsX;
        NURBS2D(2).noPtsY = noPtsY;
        NURBS2D(2).p      = p;
        NURBS2D(2).q      = q;
        NURBS2D(2).controlPts = controlPts([(noPtsX*noPtsY*noPtsZ- noPtsX*noPtsY +1) :noPtsX*noPtsY*noPtsZ],:);
        NURBS2D(2).elemsnet = elemsnet2D(noPtsX,noPtsY);
        NURBS2D(2).PID    = 2;
        NURBS2D(2).uknot  = uknot;
        NURBS2D(2).vknot  = vknot;
        NURBS2D(2).weights= weights([(noPtsX*noPtsY*noPtsZ- noPtsX*noPtsY +1) :noPtsX*noPtsY*noPtsZ]');
    
        % Face 3 : face eta = 0
        NURBS2D(3).noPtsX = noPtsX;
        NURBS2D(3).noPtsY = noPtsZ;
        NURBS2D(3).p      = p;
        NURBS2D(3).q      = r;
        B = zeros(noPtsX*noPtsZ,3);
        for k = 1:noPtsZ
           B([1:noPtsX]' + (k-1)*noPtsX, :) = controlPts([[1:noPtsX] + (k-1)*noPtsX*noPtsY]',:);
        end
        NURBS2D(3).controlPts = B;
        NURBS2D(3).elemsnet = elemsnet2D(noPtsX,noPtsZ);
        NURBS2D(3).PID    = 3;
        NURBS2D(3).uknot  = uknot;
        NURBS2D(3).vknot  = wknot;
        w = zeros(noPtsX*noPtsZ,1);
        for k = 1:noPtsZ
            w([1:noPtsX]' + (k-1)*noPtsX , :) = weights([[1:noPtsX] + (k-1)*noPtsX*noPtsY]',:);
        end
        NURBS2D(3).weights= w;
    
        % Face 4 : couter face eta = 1
        NURBS2D(4).noPtsX = noPtsX;
        NURBS2D(4).noPtsY = noPtsZ;
        NURBS2D(4).p      = p;
        NURBS2D(4).q      = r;
        B = zeros(noPtsX*noPtsZ,3);
        for k = 1:noPtsZ
           B([1:noPtsX]' + (k-1)*noPtsX, :) = controlPts([[noPtsX*noPtsY:-1:(noPtsX*noPtsY - noPtsX +1)] + (k-1)*noPtsX*noPtsY]',:);
        end
        NURBS2D(4).controlPts = B;
        NURBS2D(4).elemsnet = elemsnet2D(noPtsX,noPtsZ);
        NURBS2D(4).PID    = 4;
        NURBS2D(4).uknot  = sort(1 - uknot); 
        NURBS2D(4).vknot  = wknot;
        w = zeros(noPtsX*noPtsZ,1);
        for k = 1:noPtsZ
            w([1:noPtsX]' + (k-1)*noPtsX , :) = weights([[noPtsX*noPtsY:-1:(noPtsX*noPtsY - noPtsX +1)] + (k-1)*noPtsX*noPtsY]',:);
        end
        NURBS2D(4).weights= w;
    
        % Face 5 : left xsi = 0
        NURBS2D(5).noPtsX = noPtsY;
        NURBS2D(5).noPtsY = noPtsZ;
        NURBS2D(5).p      = q;
        NURBS2D(5).q      = r;
        B = zeros(noPtsY*noPtsZ,3);
        for k = 1:noPtsZ
           B([1:noPtsY]' + (k-1)*noPtsY, :) = controlPts([[noPtsX*noPtsY-noPtsX+1:-noPtsX:1] + (k-1)*noPtsX*noPtsY]',:);
        end
        NURBS2D(5).controlPts = B;
        NURBS2D(5).elemsnet = elemsnet2D(noPtsY,noPtsZ);
        NURBS2D(5).PID    = 5;
        NURBS2D(5).uknot  = sort(1-vknot); 
        NURBS2D(5).vknot  = wknot;
        w = zeros(noPtsY*noPtsZ,1);
        for k = 1:noPtsZ
            w([1:noPtsY]' + (k-1)*noPtsY , :) = weights([[noPtsX*noPtsY-noPtsX+1:-noPtsX:1] + (k-1)*noPtsX*noPtsY]',:);
        end
        NURBS2D(5).weights= w;
    
        % Face 6 : right xsi = 1
        NURBS2D(6).noPtsX = noPtsY;
        NURBS2D(6).noPtsY = noPtsZ;
        NURBS2D(6).p      = q;
        NURBS2D(6).q      = r;
        B = zeros(noPtsY*noPtsZ,3);
        for k = 1:noPtsZ
           B([1:noPtsY]' + (k-1)*noPtsY, :) = controlPts([[noPtsX:noPtsX:noPtsX*noPtsY] + + (k-1)*noPtsX*noPtsY]',:);
        end
        NURBS2D(6).controlPts = B;
        NURBS2D(6).elemsnet = elemsnet2D(noPtsY,noPtsZ);
        NURBS2D(6).PID    = 6;
        NURBS2D(6).uknot  = vknot; 
        NURBS2D(6).vknot  = wknot;
        w = zeros(noPtsY*noPtsZ,1);
        for k = 1:noPtsZ
            w([1:noPtsY]' + (k-1)*noPtsY , :) = weights([[noPtsX:noPtsX:noPtsX*noPtsY] + + (k-1)*noPtsX*noPtsY]',:);
        end
        NURBS2D(6).weights= w;
        
        % Save NURBS 2D as boundary of 3D NURBS
        NURBS(i).NURBS2D = NURBS2D;
    end
    NURBSwith2D=NURBS;
end
