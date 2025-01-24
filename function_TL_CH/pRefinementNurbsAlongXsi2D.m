function [Bnew,weights_new,uKnotnew,pb,noPtsX_new,T] = pRefinementNurbsAlongXsi2D(B,uKnot,noPtsX,noPtsY,p,ele_ord,weights)
    Bxsi = B([1:noPtsX],:);
    weightsXsi = weights([1:noPtsX],:);
    if length(B(1,:))==3
        a=3;
    else 
        a=2;
    end
    [~,uKnotnew,~,Txsi,pb] = pRefinement(uKnot,Bxsi,p,ele_ord,weightsXsi);
    noPtsX_new = length(uKnotnew)-pb-1;
    Bnew = zeros(noPtsX_new*noPtsY,a);
    weights_new = zeros(noPtsX_new*noPtsY,1);
    T    = zeros(noPtsX_new*noPtsY,noPtsX*noPtsY);
    for i = 1:noPtsY
        numline = (i-1)*noPtsX+1: noPtsX*i;
        numline_new =  (i-1)*noPtsX_new+1:noPtsX_new*i;
        T(numline_new,numline)=Txsi;
    end
    weights_new = T*weights;
    Bnew = T*(B.*weights);
    Bnew = Bnew ./ weights_new;
end

