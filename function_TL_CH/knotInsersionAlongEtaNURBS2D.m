function [Bnew,vKnot_new,weights_new,noPtsY_new,T]=knotInsersionAlongEtaNURBS2D(vKnot,eta_new,B,q,noPtsX,noPtsY,weights)
    noPtsY_new = noPtsY+1;
    B = B.*weights;
    Beta = B([1:noPtsX:((noPtsX*noPtsY)-noPtsX+1)],:);
    
     if length(B(1,:))==3
        a=3;
    else 
        a=2;
     end
    
    % working in progress
    [~,Teta] = knotInserstion(vKnot,eta_new,Beta,q);
    Bnew = zeros(noPtsX*noPtsY_new,a);
    weights_new = zeros(noPtsX*noPtsY_new,1);
    T    = zeros(noPtsX*noPtsY_new,noPtsX*noPtsY);
    for i = 1:noPtsX
        %numline = (1:noPtsX:(noPtsX-1)*noPtsY+1) + (i-1)
        %numline_new = (1:noPtsX:(noPtsX-1)*(noPtsY_new+1)) + (i-1)
        numline = (1:noPtsX:((noPtsX*noPtsY)-noPtsX+1)) + (i-1);
        numline_new = (1:noPtsX:((noPtsX*noPtsY_new)-noPtsX+1)) + (i-1);
        weights_new(numline_new) =  Teta*weights(numline);
        Bnew(numline_new,:) = Teta*B(numline,:);
        T(numline_new,numline)=Teta;
     end
    Bnew = Bnew./weights_new;
    [vKnot_new] = newKnot(vKnot,eta_new);
end

