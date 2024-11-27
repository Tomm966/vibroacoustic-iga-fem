function [IGA_c] = AssembleCouplingIGA(INTERFACE,IGA_s,IGA_f)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% INTERFACE: This input contains information about the interface between two domains (e.g., structure and fluid). It includes:
% INTERFACE.NURBS2D.controlPts: Control points defining the NURBS surface at the interface.
% INTERFACE.NURBS2D.Cfs: Coupling matrix at the interface between the structure and fluid.
% INTERFACE.NURBS2D.idF: Local degrees of freedom (DOF) on the fluid side.
% IGA_s: This is the IGA (Isogeometric Analysis) structure for the structural domain, which includes:
% IGA_s.B: Matrix of control points or degrees of freedom (DOFs) for the structure.
% IGA_s.nbB: Number of basis functions or control points for the structure.
% IGA_f: This is the IGA structure for the fluid domain, which includes:
% IGA_f.B: Matrix of control points or DOFs for the fluid.
% IGA_f.nbB: Number of basis functions or control points for the fluid.
% OUTPUT:
% IGA_c: This is the coupled IGA structure that contains:
% IGA_c.Cfs: The sparse coupling matrix between the structure and fluid. It has dimensions nbB_s x nbB_f, where nbB_s is the number of DOFs in the structure and nbB_f is the number of DOFs in the fluid.

    IGA_c = [];
    IGA_c.Cfs = sparse(IGA_s.nbB,IGA_f.nbB);
    [~,ind_in_common_structure] = ismembertol(INTERFACE.NURBS2D.controlPts,IGA_s.B,'ByRows',1e-6);
    INTERFACE.NURBS2D.global_s = nonzeros(ind_in_common_structure);

    [~,ind_in_common_fluid] = ismembertol(INTERFACE.NURBS2D.controlPts,IGA_f.B,'ByRows',1e-6);
    INTERFACE.NURBS2D.global_f = nonzeros(ind_in_common_fluid);
   
    

    dofuglobal = [                      INTERFACE.NURBS2D.global_s ; ...
                      size(IGA_s.B,1) + INTERFACE.NURBS2D.global_s ; ...
                    2*size(IGA_s.B,1) + INTERFACE.NURBS2D.global_s ];
    dofpglobal = [ INTERFACE.NURBS2D.global_f];

    dofplocal  = [INTERFACE.NURBS2D.idF];

    IGA_c.Cfs(dofuglobal,dofpglobal) = IGA_c.Cfs(dofuglobal,dofpglobal) + ...
                                       INTERFACE.NURBS2D.Cfs;
end

