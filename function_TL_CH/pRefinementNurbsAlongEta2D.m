function [Bnew,weights_new,vKnotnew,qb,noPtsY_new,T] = pRefinementNurbsAlongEta2D(B,vKnot,noPtsX,noPtsY,q,ele_ord_eta,weights)
    Beta = B([1:noPtsX:((noPtsX*noPtsY)-noPtsX+1)],:);
    weightseta = weights([1:noPtsX:((noPtsX*noPtsY)-noPtsX+1)],:);
    if length(B(1,:))==3
        a=3;
    else 
        a=2;
    end
    [~,vKnotnew,~,Teta,qb] = pRefinement(vKnot,Beta,q,ele_ord_eta,weightseta);
    noPtsY_new = length(vKnotnew)-qb-1;
    Bnew = zeros(noPtsX*noPtsY_new,a);
    weights_new = zeros(noPtsX*noPtsY_new,1);
    T    = zeros(noPtsX*noPtsY_new,noPtsX*noPtsY);
    for i = 1:noPtsX
        numline = (1:noPtsX:(noPtsX-1)*noPtsY+1) + (i-1);
        numline_new = (1:noPtsX:(noPtsX-1)*(noPtsY_new+1)) + (i-1);
        numline = (1:noPtsX:((noPtsX*noPtsY)-noPtsX+1)) + (i-1);
        numline_new = (1:noPtsX:((noPtsX*noPtsY_new)-noPtsX+1)) + (i-1);
        T(numline_new,numline)=Teta;
    end
    weights_new = T*weights;
    Bnew = T*(B.*weights);
    Bnew = Bnew ./ weights_new;
end

