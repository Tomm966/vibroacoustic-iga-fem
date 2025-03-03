function [NURBSnew_struct,Tnew]=refine3D(NURBS,number_of_new_xsi,...
         number_of_new_ni,number_of_new_eta)

controlPts = NURBS.controlPts;
noPtsX = NURBS.noPtsX;
noPtsY = NURBS.noPtsY;
noPtsZ = NURBS.noPtsZ;
weights = NURBS.weights;
uKnot = NURBS.uknot;
vKnot = NURBS.vknot;
wKnot = NURBS.wknot;
p = NURBS.p;
q = NURBS.q;
r = NURBS.r;
[xsi_vec_new] = equispaceKnot(uKnot,number_of_new_xsi);
% xsi_vec_new = linspace(1/(number_of_new_xsi+1),1-1/(number_of_new_xsi+1),number_of_new_xsi);
Bold = controlPts;
uKnot_old = uKnot;
noPtsX_old = noPtsX;
weights_old = weights;
for i = 1:length(xsi_vec_new)
    xsi_new = xsi_vec_new(i);
    [Bnew,uKnot_new,weights_new,noPtsX_new,T]=knotInsersionAlongXsiNURBS3D(uKnot_old,xsi_new,Bold,p,noPtsX_old,noPtsY,noPtsZ,weights_old);
    if i == 1;
        Tnew = T;
    else
        Tnew = T*Tnew;
    end
    Bold = Bnew;
    uKnot_old = uKnot_new;
    noPtsX_old = noPtsX_new;
    xsi_new = xsi_vec_new(i);
    weights_old = weights_new;
end
[ni_vec_new] = equispaceKnot(vKnot,number_of_new_ni);
% ni_vec_new = linspace(1/(number_of_new_ni+1),1-1/(number_of_new_ni+1),number_of_new_ni);
Bold = Bnew;
vKnot_old = vKnot;
noPtsY_old = noPtsY;
weights_old = weights_new;
for i = 1:length(ni_vec_new)
    ni_new = ni_vec_new(i);
    [Bnew,vKnot_new,weights_new,noPtsY_new,T]=knotInsersionAlongNiNURBS3D(vKnot_old, ...
                                                                   ni_new,...
                                                                   Bold,...
                                                                   q,...
                                                                   noPtsX_new,...
                                                                   noPtsY_old,...
                                                                   noPtsZ,weights_old);
    Bold = Bnew;
    vKnot_old = vKnot_new;
    noPtsY_old = noPtsY_new;
    ni_new = ni_vec_new(i);
    Tnew = T*Tnew;
    weights_old=weights_new;
end
[eta_vec_new] = equispaceKnot(wKnot,number_of_new_eta);
% eta_vec_new = linspace(1/(number_of_new_eta+1),1-1/(number_of_new_eta+1),number_of_new_eta);
Bold = Bnew;
wKnot_old = wKnot;
noPtsZ_old = noPtsZ;
weights_old = weights_new;
for i = 1:length(eta_vec_new)
    eta_new = eta_vec_new(i);
    [Bnew,wKnot_new,weights_new,noPtsZ_new,T]=knotInsersionAlongEtaNURBS3D(wKnot_old,eta_new,Bold,r,noPtsX_new,...
                                                                            noPtsY_new,noPtsZ_old,weights_old);
    Bold = Bnew;
    wKnot_old = wKnot_new;
    noPtsZ_old = noPtsZ_new;
    eta_new = eta_vec_new(i);
    Tnew = T*Tnew;
    weights_old = weights_new;
end

NURBSnew_struct.controlPts = Bnew;
NURBSnew_struct.noPtsX  = noPtsX_new ;
NURBSnew_struct.noPtsY  = noPtsY_new ;
NURBSnew_struct.noPtsZ  = noPtsZ_new ;
NURBSnew_struct.weights = weights_new ;
NURBSnew_struct.uknot = uKnot_new ;
NURBSnew_struct.vknot = vKnot_new ;
NURBSnew_struct.wknot = wKnot_new ;
NURBSnew_struct.p = p;
NURBSnew_struct.q = q;
NURBSnew_struct.r = r;
end
