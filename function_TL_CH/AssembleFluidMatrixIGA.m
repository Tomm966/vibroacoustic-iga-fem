function [IGA_f,NURBSnew_fluid] = AssembleFluidMatrixIGA(NURBSnew_fluid)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% NURBSnew_fluid: A structure array representing different NURBS patches for the fluid domain. Each element i of this structure contains:
% NURBSnew_fluid(i).controlPts: Control points for the i-th NURBS patch.
% NURBSnew_fluid(i).K: The local stiffness matrix for the i-th NURBS patch.
% NURBSnew_fluid(i).M: The local mass matrix for the i-th NURBS patch.
% NURBSnew_fluid(i).global_to_local: A field that will store the indices mapping the global control points to the local patch control points after the function processes them.
% OUTUP:
% IGA_f: This is the assembled fluid IGA structure which contains:
% IGA_f.Bconcatenated: A matrix concatenating all control points from the different NURBS patches.
% IGA_f.B: A matrix containing the unique control points from Bconcatenated, with duplicates removed.
% IGA_f.nbB: The number of unique control points (degrees of freedom) in the fluid domain.
% IGA_f.K: The global stiffness matrix for the fluid, constructed by assembling the local stiffness matrices from each patch.
% IGA_f.M: The global mass matrix for the fluid, constructed by assembling the local mass matrices from each patch.
% NURBSnew_fluid: This is the updated input structure, where each NURBSnew_fluid(i).global_to_local is now filled with the global indices corresponding to the local control points of the i-th patch.

%======================================
% Structure for the global system
%======================================
IGA_f = [];
IGA_f.Bconcatenated = vertcat(NURBSnew_fluid.controlPts);     % concatenate all the nodes
IGA_f.B = uniquetol(IGA_f.Bconcatenated,'ByRows',true); % find the unique nodes
for i = 1:length(NURBSnew_fluid)
    [~,ind_in_common] = ismembertol(NURBSnew_fluid(i).controlPts,IGA_f.B,'ByRows',1e-6);
    NURBSnew_fluid(i).global_to_local = nonzeros(ind_in_common);
end
%---------------------------------------
IGA_f.nbB = size(IGA_f.B,1);
IGA_f.K = sparse(IGA_f.nbB,IGA_f.nbB);
IGA_f.M = sparse(IGA_f.nbB,IGA_f.nbB);
for i = 1:length(NURBSnew_fluid)
    dofpatch = NURBSnew_fluid(i).global_to_local;
    IGA_f.K(dofpatch,dofpatch) = IGA_f.K(dofpatch,dofpatch) + NURBSnew_fluid(i).K;
    IGA_f.M(dofpatch,dofpatch) = IGA_f.M(dofpatch,dofpatch) + NURBSnew_fluid(i).M;
end
end

