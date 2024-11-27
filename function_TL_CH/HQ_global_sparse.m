function [H,Q] = HQ_global_sparse(fmesh,c0,rhof)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT
% fmesh: A structure containing mesh information, which includes:
% ndof: Total number of degrees of freedom (DOFs).
% nelem: Total number of elements in the mesh.
% elems: Element connectivity matrix, where each row corresponds to an element and contains the indices of its nodes.
% nodes: Coordinates of the nodes in the mesh, where each row represents a node's coordinates (x, y, z).
% c0: A constant value used in the local computation of H and Q. Its specific significance depends on the context of the problem.
% 
% rhof: Density or material property, again context-dependent, used in the local computation.

% OUTPUT:
% H: A sparse matrix representing the global stiffness matrix or similar (depending on context) with dimensions based on the total number of degrees of freedom.
% Q: A sparse matrix representing the global mass matrix or similar, also with dimensions based on the total number of degrees of freedom.

fndof  = fmesh.ndof  ;
fnelem = fmesh.nelem ;
felems = fmesh.elems ;
fnodes = fmesh.nodes ;


FEM = [];

for i = 1:fnelem
 
    % Global-local numerotation P --> 20 ddl per element
    %---------------------------------------------------------------------%    
    idnod = felems(i,:);
    dofp  = idnod;
    dofelems = dofp;

    %---------find(ismember(selems,ddlu));------------------------------------------------------------%
    
    X = fnodes(idnod,1);
    Y = fnodes(idnod,2);
    Z = fnodes(idnod,3);
    
    [Hl,Ql] = HQ_local(X,Y,Z,c0,rhof);
    
    state = ['H and Q asssembly ',num2str(i),' sur ',num2str(fnelem)];

    FEM(i).Ig = kron(ones(20,1),dofelems');
    FEM(i).Jg = kron(dofelems',ones(20,1));
    FEM(i).Hg = reshape(Hl,20*20,1);
    FEM(i).Qg = reshape(Ql,20*20,1);  
end

Hg = vertcat(FEM.Hg);
Qg = vertcat(FEM.Qg);
Ig = vertcat(FEM.Ig);
Jg = vertcat(FEM.Jg);
        
H = sparse(Ig,Jg,Hg,fndof,fndof);
Q = sparse(Ig,Jg,Qg,fndof,fndof);


end

