function [K,M] = KM_global_sparse(smesh,Cs,rho)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% smesh: A structure containing the mesh information for the finite element analysis, which includes:
% ndof: The total number of degrees of freedom in the mesh.
% nelem: The total number of elements in the mesh.
% elems: A matrix where each row represents an element and contains the indices of the nodes that form that element.
% nodes: A matrix where each row represents a node and contains its coordinates (x, y, z).
% Cs: A constant or matrix related to the material properties or stiffness coefficients used in the computation of the global stiffness matrix.
% rho: A scalar representing the material density, used in the computation of the global mass matrix.

% Output:
% K: A sparse matrix representing the global stiffness matrix, of size (sndof x sndof).
% M: A sparse matrix representing the global mass matrix, also of size (sndof x sndof).

%K = sparse(zeros(sndof));
%M = sparse(zeros(sndof));

sndof  = smesh.ndof  ;
snelem = smesh.nelem ;
selems = smesh.elems ;
snodes = smesh.nodes ;

FEM = [];

%--- Initialize Gauss pooints -----------------------------------------
    
intPt = [];
[ir wp N dNeta dNnu dNte] = Point_intreg_3D(); 
intPt.ir = ir;    
intPt.wp = wp;
intPt.N  = N ;
intPt.dNeta = dNeta ;
intPt.dNnu  = dNnu  ;
intPt.dNte  = dNte  ;

for i = 1:snelem
 
    % Global-local numerotation U --> 60 ddl per element
    %---------------------------------------------------------------------%    
    idnod = selems(i,:);
    dofx  = 3*idnod-2;
    dofy  = 3*idnod-1;
    dofz  = 3*idnod;
    dofelems_concat = [dofx' dofy' dofz'];
    dofelems = reshape(dofelems_concat',[1,60]);
    %---------------------------------------------------------------------%
    
    X = snodes(idnod,1);
    Y = snodes(idnod,2);
    Z = snodes(idnod,3);

    [Kl Ml] = KM_local(X,Y,Z,Cs,rho,intPt);

    FEM(i).Ig = kron(ones(60,1),dofelems');
    FEM(i).Jg = kron(dofelems',ones(60,1));
    FEM(i).Kg = reshape(Kl,60*60,1);
    FEM(i).Mg = reshape(Ml,60*60,1);        


end

Kg = vertcat(FEM.Kg);
Mg = vertcat(FEM.Mg);
Ig = vertcat(FEM.Ig);
Jg = vertcat(FEM.Jg);
        
K = sparse(Ig,Jg,Kg,sndof,sndof);
M = sparse(Ig,Jg,Mg,sndof,sndof);


end


