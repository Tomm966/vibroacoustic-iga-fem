function [IGA_s,NURBSnew_struct] = AssembleStructMatrixIGA(NURBSnew_struct)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% NURBSnew_struct: A structure array that represents different NURBS patches for the structural domain. Each element i in this structure contains:
% NURBSnew_struct(i).controlPts: Control points of the i-th NURBS patch for the structure.
% NURBSnew_struct(i).K: The local stiffness matrix for the i-th NURBS patch.
% NURBSnew_struct(i).M: The local mass matrix for the i-th NURBS patch.
% NURBSnew_struct(i).global_to_local: A field that will store the indices mapping the global control points to the local control points after processing.
% OUTPUT:
% IGA_s: This is the assembled IGA structure for the structural domain, which contains:
% IGA_s.Bconcatenated: A matrix of concatenated control points from all NURBS patches for the structure.
% IGA_s.B: A matrix containing the unique control points from Bconcatenated, with duplicates removed.
% IGA_s.nbB: The number of degrees of freedom (DOFs) in the structure domain. Since there are 3 DOFs per control point (likely representing x, y, and z directions), this is set to 3 * size(IGA_s.B,1).
% IGA_s.K: The global stiffness matrix for the structural domain, constructed by assembling the local stiffness matrices from each patch.
% IGA_s.M: The global mass matrix for the structural domain, assembled similarly from the local mass matrices.
% NURBSnew_struct: This is the updated input structure where the field NURBSnew_struct(i).global_to_local is now filled with the global indices corresponding to the local control points of the i-th patch.

%======================================
% Structure for the global system
%======================================
IGA_s = [];
IGA_s.Bconcatenated = vertcat(NURBSnew_struct.controlPts);     % concatenate all the nodes
IGA_s.B = uniquetol(IGA_s.Bconcatenated,'ByRows',true); % find the unique nodes
%---------------------------------------
IGA_s = [];
IGA_s.Bconcatenated = vertcat(NURBSnew_struct.controlPts);     % concatenate all the nodes
IGA_s.B = uniquetol(IGA_s.Bconcatenated,'ByRows',true); % find the unique nodes
for i = 1:length(NURBSnew_struct)
    [~,ind_in_common] = ismembertol(NURBSnew_struct(i).controlPts,IGA_s.B,'ByRows',1e-6);
    NURBSnew_struct(i).global_to_local = nonzeros(ind_in_common);
end
IGA_s.nbB = 3*size(IGA_s.B,1);
IGA_s.K = sparse(IGA_s.nbB,IGA_s.nbB);
IGA_s.M = sparse(IGA_s.nbB,IGA_s.nbB);

for i=1:length(NURBSnew_struct)

    dofpatch = [NURBSnew_struct(i).global_to_local; ...
                NURBSnew_struct(i).global_to_local+size(IGA_s.B,1); ...
                NURBSnew_struct(i).global_to_local+2*size(IGA_s.B,1)];
    
    IGA_s.K(dofpatch,dofpatch) = IGA_s.K(dofpatch,dofpatch) + NURBSnew_struct(i).K;
    IGA_s.M(dofpatch,dofpatch) = IGA_s.M(dofpatch,dofpatch) + NURBSnew_struct(i).M;
end
doftot_s=[1:IGA_s.nbB];
end

