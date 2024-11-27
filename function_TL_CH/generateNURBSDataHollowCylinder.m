function [NURBS] = generateNURBSDataHollowCylinder(R,t,L)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% R: Outer radius of the hollow cylinder.
% t: Thickness of the cylinder wall.
% L: Length of the cylinder.

% OUTPUT:
% NURBS: A structure containing data for four NURBS solids. Each solid includes:
% uknot: Knot vector in the u direction.
% vknot: Knot vector in the v direction.
% wknot: Knot vector in the w direction.
% noPtsX: Number of control points in the u direction.
% noPtsY: Number of control points in the v direction.
% noPtsZ: Number of control points in the w direction.
% weights: Weights for the control points.
% controlPts: Control points for the solid.

    Ri     = R -t;
    Rii    = 0.5*(R+Ri);
    hL     = 0.5*L;
    p = 2; %polynomial degree
    q = 2; %polynomial degree
    r = 2; %polynomial degree
    NURBS = [];
    NURBS.p = p;
    NURBS.q = q;
    NURBS.r = r;
    %----- First volume --------------------------------------------------%
    uKnot_1 = [0 0 0 1 1 1]; %starting knot vector %%xsi
    vKnot_1 = [0 0 0 1 1 1]; %starting knot vector %%ni
    wKnot_1 = [0 0 0 1 1 1]; %starting knot vector %%eta
    noPtsX_1      = length(uKnot_1)-p-1;
    noPtsY_1      = length(vKnot_1)-q-1;
    noPtsZ_1      = length(wKnot_1)-r-1;
    weights_1 = ones(1,noPtsX_1*noPtsY_1*noPtsZ_1)';
    fac = 1/sqrt(2);
    weights_1([2,5,8,11,14,17,20,23,26]) = fac;
    B_1= [ -Ri  0   0;-Ri  Ri  0;0  Ri  0;
           -Rii 0   0;-Rii Rii 0;0  Rii 0;
           -R   0   0;-R   R   0;0   R   0;
           ...
           -Ri  0  hL;-Ri  Ri  hL; 0  Ri  hL;
           -Rii 0   hL;-Rii Rii hL; 0  Rii hL;
           -R   0   hL;-R   R   hL; 0   R   hL;
           ...
           -Ri  0   L;-Ri  Ri  L; 0  Ri  L;
           -Rii 0   L;-Rii Rii L; 0  Rii L;
           -R   0   L;-R   R   L;0   R   L];
    %---- Second volume --------------------------------------------------%
    uKnot_2 = [0 0 0 1 1 1]; %starting knot vector %%xsi
    vKnot_2 = [0 0 0 1 1 1]; %starting knot vector %%ni
    wKnot_2 = [0 0 0 1 1 1]; %starting knot vector %%eta
    noPtsX_2      = length(uKnot_2)-p-1;
    noPtsY_2      = length(vKnot_2)-q-1;
    noPtsZ_2      = length(wKnot_2)-r-1;
    weights_2 = ones(1,noPtsX_2*noPtsY_2*noPtsZ_2)';
    fac = 1/sqrt(2);
    weights_2([2,5,8,11,14,17,20,23,26]) = fac;
    B_2=[0   Ri 0 ; Ri  Ri  0;Ri  0 0;
         0   Rii 0; Rii Rii 0;Rii 0 0;
         0   R 0  ; R   R   0;R   0 0;
            ...
            0   Ri hL; Ri Ri hL;Ri 0 hL;
            0  Rii hL;Rii Rii hL;Rii 0 hL;
            0  R hL;R  R hL;R   0 hL;
            ...
             0   Ri L; Ri Ri L;Ri 0 L;
            0  Rii L;Rii Rii L;Rii 0 L;
            0  R L;R  R L;R   0 L;
           ]; 
    %---- Third volume  --------------------------------------------------%
    uKnot_3 = [0 0 0 1 1 1]; %starting knot vector %%xsi
    vKnot_3 = [0 0 0 1 1 1]; %starting knot vector %%ni
    wKnot_3 = [0 0 0 1 1 1]; %starting knot vector %%eta
    p = 2; %polynomial degree
    q = 2; %polynomial degree
    r = 2; %polynomial degree
    noPtsX_3      = length(uKnot_3)-p-1;
    noPtsY_3      = length(vKnot_3)-q-1;
    noPtsZ_3      = length(wKnot_3)-r-1;
    weights_3 = ones(1,noPtsX_3*noPtsY_3*noPtsZ_3)';
    fac = 1/sqrt(2);
    weights_3([2,5,8,11,14,17,20,23,26]) = fac;
    B_3=[0   -Ri 0; -Ri -Ri 0;-Ri    0 0;
        0  -Rii 0;-Rii -Rii 0;-Rii 0 0;
        0  -R 0;-R  -R 0;-R   0 0;
        ...
        0   -Ri hL; -Ri -Ri hL;-Ri 0 hL;
        0  -Rii hL;-Rii -Rii hL;-Rii 0 hL;
        0  -R hL;-R  -R hL;-R   0 hL;
        ...
         0   -Ri L; -Ri -Ri L;-Ri 0 L;
        0  -Rii L;-Rii -Rii L;-Rii 0 L;
        0  -R L;-R  -R L;-R   0 L;
       ];
    %---- Fourth volume --------------------------------------------------%
    uKnot_4 = [0 0 0 1 1 1]; %starting knot vector %%xsi
    vKnot_4 = [0 0 0 1 1 1]; %starting knot vector %%ni
    wKnot_4 = [0 0 0 1 1 1]; %starting knot vector %%eta
    p = 2; %polynomial degree
    q = 2; %polynomial degree
    r = 2; %polynomial degree
    noPtsX_4      = length(uKnot_4)-p-1;
    noPtsY_4      = length(vKnot_4)-q-1;
    noPtsZ_4      = length(wKnot_4)-r-1;
    weights_4 = ones(1,noPtsX_4*noPtsY_4*noPtsZ_4)';
    fac = 1/sqrt(2);
    weights_4([2,5,8,11,14,17,20,23,26]) = fac;
    B_4= [ Ri  0   0;Ri  -Ri  0;0  -Ri  0;
             Rii 0   0;Rii -Rii 0;0  -Rii 0;
             R   0   0;R   -R   0;0   -R   0;
              ...
              Ri  0  hL;Ri  -Ri  hL; 0  -Ri  hL;
             Rii 0   hL;Rii -Rii hL; 0  -Rii hL;
             R   0   hL;R   -R   hL; 0   -R   hL;
             ...
             Ri  0   L;Ri  -Ri  L; 0  -Ri  L;
             Rii 0   L;Rii -Rii L; 0  -Rii L;
             R   0   L;R   -R   L;0   -R   L];
    %---- Solid 1 --------------------------------------------------------%
    NURBS(1).uknot = uKnot_1;
    NURBS(1).vknot = vKnot_1;
    NURBS(1).wknot = wKnot_1;
    NURBS(1).noPtsX = noPtsX_1;
    NURBS(1).noPtsY = noPtsY_1;
    NURBS(1).noPtsZ = noPtsZ_1;
    NURBS(1).weights= weights_1;
    NURBS(1).controlPts = B_1;
    NURBS(1).p = p;
    NURBS(1).q = q;
    NURBS(1).r = r;

    %---- Solid 2 --------------------------------------------------------%
    NURBS(2).uknot = uKnot_2;
    NURBS(2).vknot = vKnot_2;
    NURBS(2).wknot = wKnot_2;
    NURBS(2).noPtsX = noPtsX_2;
    NURBS(2).noPtsY = noPtsY_2;
    NURBS(2).noPtsZ = noPtsZ_2;
    NURBS(2).weights= weights_2;
    NURBS(2).controlPts = B_2;
    NURBS(2).p = p;
    NURBS(2).q = q;
    NURBS(2).r = r;

    %---- Solid 3 --------------------------------------------------------%
    NURBS(3).uknot = uKnot_3;
    NURBS(3).vknot = vKnot_3;
    NURBS(3).wknot = wKnot_3;
    NURBS(3).noPtsX = noPtsX_3;
    NURBS(3).noPtsY = noPtsY_3;
    NURBS(3).noPtsZ = noPtsZ_3;
    NURBS(3).weights= weights_3;
    NURBS(3).controlPts = B_3;
    NURBS(3).p = p;
    NURBS(3).q = q;
    NURBS(3).r = r;

    %---- Solid 4 --------------------------------------------------------%
    NURBS(4).uknot = uKnot_4;
    NURBS(4).vknot = vKnot_4;
    NURBS(4).wknot = wKnot_4;
    NURBS(4).noPtsX = noPtsX_4;
    NURBS(4).noPtsY = noPtsY_4;
    NURBS(4).noPtsZ = noPtsZ_4;
    NURBS(4).weights= weights_4;
    NURBS(4).controlPts = B_4;
    NURBS(4).p = p;
    NURBS(4).q = q;
    NURBS(4).r = r;
end

