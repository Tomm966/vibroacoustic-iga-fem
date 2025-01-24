function [Bnew,uKnot_new,weights_new,noPtsX_new,T] = knotInsersionAlongXsiNURBS2D(uKnot,xsi_new,B,p,noPtsX,noPtsY,weights)
    noPtsX_new = noPtsX+1;
    B = B.*weights;
    Bxsi = B([1:noPtsX],:);
    
    if length(B(1,:))==3
        a=3;
    else 
        a=2;
    end
    [~,Txsi] = knotInserstion(uKnot,xsi_new,Bxsi,p);
    Bnew = zeros(noPtsX_new*noPtsY,a);
    weights_new = zeros(noPtsX_new*noPtsY,1);
    T    = zeros(noPtsX_new*noPtsY,noPtsX*noPtsY);
    for i = 1:noPtsY
        numline = (i-1)*noPtsX+1: noPtsX*i;
        numline_new =  (i-1)*noPtsX_new+1:noPtsX_new*i;
        weights_new(numline_new) =  Txsi*weights(numline);
        Bnew(numline_new,:) = Txsi*B(numline,:);
        T(numline_new,numline)=Txsi;
    end
    Bnew = Bnew./weights_new;
    [uKnot_new] = newKnot(uKnot,xsi_new);
end

