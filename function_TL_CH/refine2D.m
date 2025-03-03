function [NURBSnew,T]=refine2D(NURBS2D,number_of_new_xsi,number_of_new_eta)

B = NURBS2D.controlPts;
noPtsX = NURBS2D.noPtsX;
noPtsY = NURBS2D.noPtsY;
weights = NURBS2D.weights;
uKnot = NURBS2D.uknot;
vKnot = NURBS2D.vknot;
p = NURBS2D.p;
q = NURBS2D.q;

xsi_vec_new = equispaceKnot(uKnot,number_of_new_xsi);
% xsi_vec_new = linspace(1/(number_of_new_xsi+1),1-1/(number_of_new_xsi+1),number_of_new_xsi);
Bold = B;
uKnot_old = uKnot;
noPtsX_old = noPtsX;
weights_old= weights;
for i = 1:length(xsi_vec_new)
    xsi_new = xsi_vec_new(i);
    [Bnew,uKnot_new,weights_new,noPtsX_new,T]=knotInsersionAlongXsiNURBS2D(uKnot_old,xsi_new,Bold,p,noPtsX_old,noPtsY,weights_old);
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
eta_vec_new = equispaceKnot(vKnot,number_of_new_eta);
% eta_vec_new = linspace(1/(number_of_new_eta+1),1-1/(number_of_new_eta+1),number_of_new_eta);
Bold = Bnew;
vKnot_old = vKnot;
noPtsY_old = noPtsY;
weights_old = weights_new;
for i = 1:length(eta_vec_new)
    eta_new = eta_vec_new(i);
    [Bnew,vKnot_new,weights_new,noPtsY_new,T]=knotInsersionAlongEtaNURBS2D(vKnot_old,eta_new,Bold,q,noPtsX_new,noPtsY_old,weights_old);
    Bold = Bnew;
    vKnot_old = vKnot_new;
    noPtsY_old = noPtsY_new;
    eta_new = eta_vec_new(i);
    Tnew = T*Tnew;
    weights_old = weights_new;
end

NURBSnew.controlPts = Bnew;
NURBSnew.noPtsX  = noPtsX_new ;
NURBSnew.noPtsY  = noPtsY_new ;
NURBSnew.weights = weights_new ;
NURBSnew.uknot = uKnot_new ;
NURBSnew.vknot = vKnot_new ;
NURBSnew.p = p;
NURBSnew.q = q;

end
