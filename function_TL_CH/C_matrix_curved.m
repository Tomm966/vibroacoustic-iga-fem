function [Cfs, PhyGP, Normals,Area,idControlPoint_Coupl_in_fluid,...
        idControlPoint_Coupl_in_struct] = C_matrix_curved(controlPts_coupl,noCtrPts_struct,noCtrPts_fluid,controlPts_struct,controlPts_fluid, ...
                          uKnot_coupl,vKnot_coupl,noPtsX_coupl,noPtsY_coupl,p,q,weights_coupl,noGPs)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% controlPts_coupl: Control points along the coupling interface.
% noCtrPts_struct: Number of control points in the structural domain.
% noCtrPts_fluid: Number of control points in the fluid domain.
% controlPts_struct: Control points in the structural domain.
% controlPts_fluid: Control points in the fluid domain.
% uKnot_coupl: Knot vector in the u-direction for the coupling surface.
% vKnot_coupl: Knot vector in the v-direction for the coupling surface.
% noPtsX_coupl: Number of control points in the u-direction on the coupling surface.
% noPtsY_coupl: Number of control points in the v-direction on the coupling surface.
% p: Degree of NURBS basis functions in the u-direction.
% q: Degree of NURBS basis functions in the v-direction.
% weights_coupl: Weights associated with the NURBS basis functions.
% noGPs: Number of Gauss points for integration.
% OUTPUT:
% Cfs: Coupling matrix between the structural and fluid domains.
% PhyGP: Physical coordinates of Gauss points on the interface.
% Normals: Normal vectors at the Gauss points.
% Area: Total area of the coupling interface.
% idControlPoint_Coupl_in_fluid: Indices of the interface control points in the fluid domain.
% idControlPoint_Coupl_in_struct: Indices of the interface control points in the structural domain.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Elements coupling matrix
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[element_coupl,elRangeU_coupl,...
 elConnU_coupl,elRangeV_coupl,...
 elConnV_coupl,noElems_coupl,index_coupl] = generateIGA2DMesh_CH(uKnot_coupl, ...
                                        vKnot_coupl, ...
                                        noPtsX_coupl,...
                                        noPtsY_coupl, ...
                                        p,q);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Identify controle points in interface
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the coupled nodes
% nb_control_point_coupl = size(controlPts_coupl,1);
% idControlPoint_Coupl_in_fluid = zeros(nb_control_point_coupl,1);
% idControlPoint_Coupl_in_struct = zeros(nb_control_point_coupl,1);

% Iteration of fluid and structural control points


tol = 10^-8;
[~,idControlPoint_Coupl_in_fluid]=ismembertol(controlPts_coupl(:,:),controlPts_fluid(:,1:3),'Byrow',tol);
[~,idControlPoint_Coupl_in_struct]=ismembertol(controlPts_coupl(:,:),controlPts_struct(:,1:3),'Byrows',tol);


Cfs = sparse(3*noCtrPts_struct,noCtrPts_fluid);      % global stiffness matrix structure

jacob   = zeros(2,2);
Nxi     = zeros(1,p+1);
Neta    = zeros(1,q+1);
dNdxi   = zeros(1,p+1);
dNdeta  = zeros(1,q+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equations
% Stiffness matrix and external force vector

% disp('  ASSEMBLING THE SYSTEM')

% Loop over elements (knot spans)

Normals = [];
PhyGP   = [];
Area    = 0;

for e=1:noElems_coupl


   idu    = index_coupl(e,1);
   idv    = index_coupl(e,2);
   xiE    = elRangeU_coupl(idu,:); % [xi_i,xi_i+1]
   etaE   = elRangeV_coupl(idv,:); % [eta_j,eta_j+1]
   connU  = elConnU_coupl(idu,:);
   connV  = elConnV_coupl(idv,:);
   
   noFnsU = length(connU);
   noFnsV = length(connV);

   sctr = element_coupl(e,:);
   sctr_struct_x  = idControlPoint_Coupl_in_struct(element_coupl(e,:))';
   sctr_struct_y  = idControlPoint_Coupl_in_struct(element_coupl(e,:))' + noCtrPts_struct;
   sctr_struct_z  = idControlPoint_Coupl_in_struct(element_coupl(e,:))' + 2*noCtrPts_struct;
   
   sctr_struct = [sctr_struct_x,sctr_struct_y,sctr_struct_z];
   sctr_fluid   = idControlPoint_Coupl_in_fluid(element_coupl(e,:))';          %  element scatter vector
   nn_struct    = length(sctr_struct);
   nn_fluid     = length(sctr_fluid); 

   nn     = length(sctr_struct_x); 
   
   Ns = zeros(3,nn_struct);
   Nf = zeros(1,nn_fluid);
  
 
   % loop over Gauss points 
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                            
      
      % compute coords in parameter space
      Xi      = parent2ParametricSpace(xiE,pt(1));
      Eta     = parent2ParametricSpace(etaE,pt(2)); 
      J2      = jacobianPaPaMapping(xiE,etaE);
      
      dRdxi   = [];
      dRdeta  = [];
      Rd      = [];
        
      % compute derivative of basis functions w.r.t parameter coord
      
      for in=1:noFnsU
       [Ni,dNi]  = NURBSbasis (connU(in),p,Xi,uKnot_coupl,weights_coupl); %stanno in C_files NURBSbasis.c
       Nxi(in)    = Ni;
       dNdxi(in)  = dNi;
      end
      
      for in=1:noFnsV
       [Ni,dNi]  = NURBSbasis (connV(in),q,Eta,vKnot_coupl,weights_coupl);
       Neta(in)   = Ni;
       dNdeta(in) = dNi;
      end
      
      % derivate of R=Nxi*Neta w.r.t xi and eta
      % this is derivative of shape functions in FEM
      % the following code is for B-spline only!!!
      
      for j=1:noFnsV
          for i=1:noFnsU
              dRdxi  = [dRdxi  dNdxi(i) * Neta(j)];
              dRdeta = [dRdeta Nxi(i)   * dNdeta(j)];
              Rd     = [Rd     Nxi(i)   * Neta(j)];
          end
      end
      
      % for NURBS use the following
      
      %dNdxi 
      %dNdeta
      
      %[dRdxi dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights'); 
      
      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
            
      pts = controlPts_coupl(sctr,:);
      
      % Jacobian matrix
      
      jacob(1,1) = dRdxi  * pts(:,1);
      jacob(1,2) = dRdeta * pts(:,1);
      jacob(2,1) = dRdxi  * pts(:,2);
      jacob(2,2) = dRdeta * pts(:,2);
      jacob(3,1) = dRdxi  * pts(:,3);
      jacob(3,2) = dRdeta * pts(:,3);
      
      detJxy = jacob(1,1)*jacob(2,2) - jacob(2,1)*jacob(1,2);
      detJyz = jacob(2,1)*jacob(3,2) - jacob(3,1)*jacob(2,2);
      detJzx = jacob(3,1)*jacob(1,2) - jacob(1,1)*jacob(3,2);    
      J1 = sqrt(detJxy^2 + detJyz^2 + detJzx^2);

      %J1         = det(jacob);
      
      % Jacobian inverse and spatial derivatives
      
      dRdx       = [dRdxi' dRdeta']/jacob;

      % B matrix
      %        _                                      _
      %        |  N_1,x  N_2,x  ...      0      0  ... |
      %  B  =  |      0      0  ... N_1,y  N_2,y  ... |
      %        |  N_1,y  N_2,y  ... N_1,x  N_2,x  ... |
      %        -  
      %        |  N_1  N_2  ...      0      0    ... |
      %  N  =  |    0    0  ...      N_1    N_2  ... |
      %
      %
      %  Ba =   |  N_1,x  N_2,x  ...     |
      %         |  N_1,y  N_2,y  ...     |
      %        -  
      %  Na =   |  N_1   N_2   ... |
      %


      % normal evaluated for a cylinder
      phi_pts = Rd(:)'*pts;
      xphi = phi_pts(1);
      yphi = phi_pts(2);
      zphi = phi_pts(3);
      nphi = [cos(atan2(yphi,xphi)) ; sin(atan2(yphi,xphi))];

%       disp(["Reference normal ", num2str(n') ]);
%       disp(["Computed normal ", num2str(nphi') ]);

      Ns(1 ,1:nn)        =    Rd(:)';
      Ns(2 ,nn+1:2*nn)   =    Rd(:)';
      Ns(3 ,2*nn+1:3*nn) =    Rd(:)';
      Nf                 =    Rd(:)';
    
      % normal computed with derivatives

      x_dxi = dRdxi*pts(:,1);
      y_dxi = dRdxi*pts(:,2);
      z_dxi = dRdxi*pts(:,3);
      X_xsi = [x_dxi;y_dxi;z_dxi];
      x_deta = dRdeta*pts(:,1);
      y_deta = dRdeta*pts(:,2);
      z_deta = dRdeta*pts(:,3);
      X_eta = [x_deta;y_deta;z_deta];

      n = - cross(X_xsi,X_eta) ./ norm(cross(X_xsi,X_eta),2);

      PhyGP = [PhyGP, phi_pts'];
      Normals = [Normals, n];
      
      Area = Area + J1 * J2 * wt;
      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      
      Cfs(sctr_struct,sctr_fluid)   = Cfs(sctr_struct,sctr_fluid) + Ns'  *n*  Nf * J1 * J2 * wt;
    end

    

end
end