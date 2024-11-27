function [Cfsi,Area]=couplingMatrix3DMatching(elemsfsif,elemsfsis,ndoff,ndofs,nodess,nodesf,rhof)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% elemsfsif: An array representing the fluid elements at the fluid-solid interface. Each row contains the indices of the nodes in the fluid mesh that form an element.
% elemsfsis: An array representing the solid elements at the fluid-solid interface. Each row contains the indices of the nodes in the solid mesh that form an element.
% ndoff: The total number of degrees of freedom (DOFs) in the fluid domain.
% ndofs: The total number of degrees of freedom (DOFs) in the solid domain.
% nodess: A matrix containing the coordinates of the nodes in the solid domain.
% nodesf: A matrix containing the coordinates of the nodes in the fluid domain.
% rhof: A scalar representing the density of the fluid

% OUTPUT:
% Cfsi: A sparse coupling matrix that relates the DOFs of the solid domain to the DOFs of the fluid domain. This matrix is used to describe the interaction between the solid and fluid at the interface.
% Area: A scalar that represents the total area of the fluid-solid interface, which is computed based on the geometry of the solid elements at the interface.

nelemfsif = size(elemsfsif,1);
nelemfsis = size(elemsfsis,1);
% Cfsi = zeros(ndofs,ndoff);
Cfsi = spalloc(ndofs,ndoff,500); %when the dimension is too big
Area = 0;


for i = 1:nelemfsis
    % Global-local numerotation solid
    %---------------------------------------------------------------------% 
    % disp(i)
    sidnod = elemsfsis(i,:);
    tol = 10^-15;
    [~,fcdof]=ismembertol(nodess(sidnod,:),nodesf,'ByRows',tol); % "fcdof" stands for "fluid coupling dof"
    % keyboard

    sidnodx = 3*sidnod-2    ;
    sidnody = 3*sidnod-1    ;
    sidnodz = 3*sidnod      ;

    sidnod_concat = [sidnodx',sidnody',sidnodz'];
    scdof = reshape(sidnod_concat',[1,24]);


    dofq  = scdof;
    dofp  = fcdof;
    sdofelems = dofq';   
    fdofelems = dofp';

    %---------------------------------------------------------------------%

    X = nodess(sidnod,1);
    Y = nodess(sidnod,2);
    Z = nodess(sidnod,3); 

    [Cl,Se] = C_local_courbe(X,Y,Z,rhof);
    

    Cfsi(sdofelems,fdofelems) = Cfsi(sdofelems,fdofelems) + Cl; 
    Area = Area + Se;
end

end