function [NURBSnew] = pRefine3D(NURBS,ele_ord_p,ele_ord_q,ele_ord_r)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% NURBS: An array of NURBS surfaces. Each surface is represented as a struct containing control points, knot vectors, weights, and their respective degrees in the three dimensions.
% ele_ord_p: The order of elevation for the u-direction.
% ele_ord_q: The order of elevation for the v-direction.
% ele_ord_r: The order of elevation for the w-direction.
% OUTPUT
% NURBSnew: A new array of NURBS surfaces that have been elevated in their polynomial degrees.

     NURBSnew = [];
    for i = 1:length(NURBS)
        B = NURBS(i).controlPts;
        noPtsX = NURBS(i).noPtsX;
        noPtsY = NURBS(i).noPtsY;
        noPtsZ = NURBS(i).noPtsZ;
        weights = NURBS(i).weights;
        uKnot = NURBS(i).uknot;
        vKnot = NURBS(i).vknot;
        wKnot = NURBS(i).wknot;
        p = NURBS(i).p;
        q = NURBS(i).q;
        r = NURBS(i).r;
        tic
        [Bnew,wKnotnew,weights_new,noPtsZ_new,rb,T] = pRefinementAlongEtaNURBS3D(wKnot,B,r,ele_ord_r,noPtsX,noPtsY,noPtsZ,weights);
        [Bnew_1,uKnotnew,weights_new_1,noPtsX_new,pb,T] = pRefinementAlongXsiNURBS3D(uKnot,Bnew,p,ele_ord_p,noPtsX,noPtsY,noPtsZ_new,weights_new);
        [Bnew_2,vKnotnew,weights_new_2,noPtsY_new,qb,Tnew] = pRefinementAlongNiNURBS3D(vKnot,Bnew_1,q,ele_ord_q,noPtsX_new,noPtsY,noPtsZ_new,weights_new_1);
        timeOrderElevation = toc;
        NURBSnew(i).controlPts = Bnew_2;
        NURBSnew(i).noPtsX  = noPtsX_new ;
        NURBSnew(i).noPtsY  = noPtsY_new ;
        NURBSnew(i).noPtsZ  = noPtsZ_new ;
        NURBSnew(i).weights = weights_new_2 ;
        NURBSnew(i).uknot = uKnotnew ;
        NURBSnew(i).vknot = vKnotnew ;
        NURBSnew(i).wknot = wKnotnew ;
        NURBSnew(i).p = pb;
        NURBSnew(i).q = qb;
        NURBSnew(i).r = rb;
        NURBSnew(i).T = Tnew; 
%         NURBSnew(i).timeOrderElevation = timeOrderElevation;
    end
end
