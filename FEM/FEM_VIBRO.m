%FEM vibroacoustic 
clc
clear all
close all
%% Recall all functions
addpath ../C_files/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../function_TL_CH/
%%
%-------------STRUCTURE-------------%
% %====================================================%%
% % Generate NURBS data (hollow structural cylinder)   %%
% %====================================================%%
R_s    = 1;    %external radius
L      = 3;    %length
t      = 0.05; %thickness
Rint   = R_s-t;%internal radius
[NURBS_s] = generateNURBSDataHollowCylinder(R_s,t,L); %Generate structure
%-------------------------------
nex_s =  3 ;  %NUMBER OF ELEMENTS ALONG THE CIRCUMFERENTIAL DIRECTION
ney_s =  2 ;  %NUMBER OF ELEMENTS ALONG THE THICKNESS DIRECTION
nez_s =  3 ;  %NUMBER OF ELEMENTS ALONG THE AXIAL DIRECTION
%-------------------------------
elemtype = 2 ; % quadratic element
[nodess,elemss,elems_skins]=generateMeshFromNURBS(nex_s,ney_s,nez_s,elemtype,NURBS_s); %Generate FEM mesh
%------------PLOT FEM MESH---------%
matlab2VTK_fluid(nodess(:,[1,2,3]),elemss,nodess(:,3),'MESH/Cylinder_Structure')
smesh = [];
smesh.nodes = nodess; %nodes structure
smesh.elems = elemss; %elements structure 
smesh.ndof =  3*size(nodess,1); %dofs structure
smesh.nelem = size(elemss,1);
%%
%--------STRUCTURAL PROPERTIES-------%
E     = 7e10;  %Young Modulus
nu    = 0.3;   %Poisson
rho   = 2700;  %density
lamb  = nu*E/((1+nu)*(1-2*nu));
mu    = E/(2*(1+nu));
C_str =[lamb+2*mu lamb      lamb      0  0  0;
        lamb      lamb+2*mu lamb      0  0  0;
        lamb      lamb      lamb+2*mu 0  0  0;
        0         0         0         mu 0  0;
        0         0         0         0  mu 0;
        0         0         0         0  0  mu];
%----Genrate STRUCTURAL MATRICES
[K,M] = KM_global_sparse(smesh,C_str,rho);
%----Boundary conditions (the numerotation is x1,y1,z1,x2,y2,z2,...,xend,yend,zend,
tol = 10^-6;
z0  = 0;
zL  = L;
numclamp_0 = find(ismembertol(nodess(:,3),z0,tol));
numclamp_L = find(ismembertol(nodess(:,3),zL,tol));
dofx_0 = 3*numclamp_0 - 2;
dofy_0 = 3*numclamp_0 - 1;
dofx_L = 3*numclamp_L - 2;
dofy_L = 3*numclamp_L - 1;
dofx   = vertcat(dofx_0,dofx_L);
dofy   = vertcat(dofy_0,dofy_L);
dofk = [dofx;dofy];          %blocked dofs
doftot = [1:smesh.ndof]';
dofu = setdiff(doftot,dofk); %free dofs
%%
%----Structural Modal Analysis
nmode_solid = 150; %number of modes to search
[vecp_s,valp_s]  = eigs(K(dofu,dofu),M(dofu,dofu),nmode_solid,'smallestreal','tolerance',0.08,'SubspaceDimension',2*nmode_solid+35,'MaxIterations',500);
[freq_solid,renums] = sort(sqrt(diag(real(valp_s)))/(2*pi));
freq_solid=real(freq_solid); %natural frequencies
disp([freq_solid(1:10)]) 
Phi_s = zeros(smesh.ndof,nmode_solid); %eigenvectors
Phi_s(dofu,:) = real(vecp_s);
%-----PLOT Structural modes
for i = 1:20 %number of mode to visualize in Paraview
    matlab2VTK_structure(nodess,elemss,Phi_s(:,i),['Structural_Mode/Modestructure_',num2str(i)])  
end
%%
%-------------CAVITY-------------%
%----Generate NURBS datas
R_f      = R_s-t; %fluid radius
[NURBS_f] = generateNURBSDataCylinder(R_f,L); %generate fluid cavity
ney_f =  nex_s ;
[nodesf,elemsf,elems_skinf]=generateMeshFromNURBS(nex_s,ney_f,nez_s,2,NURBS_f);
%------PLOT FLUID FEM MESH
matlab2VTK_fluid(nodesf(:,[1,2,3]),elemsf,nodesf(:,3),'MESH/Cylinder_Cavity')
fmesh = [];
fmesh.nodes = nodesf(:,[1,2,3]); %fluid nodes
fmesh.elems = elemsf; %elements nodes
fmesh.ndof  = size(nodesf,1); %dofs nodes
fmesh.nelem = size(elemsf,1);
%%
%--------ACOUSTIC PROPERTIES-------%
c0   = 340.0; % FLUID SPEED (m/s)
rhof = 1.2;   % FLUID DENSITY (kg/m^3)
%---BUILD ACOUSTIC MATRICES
[H,Q] = HQ_global_sparse(fmesh,c0,rhof);
%---ACOUSTIC MODAL ANALYSIS
nmode_fluid = 150; %fluid mode to obtain
[vecp_f,valp_f] = eigs(H,Q,nmode_fluid,'smallestreal','tolerance',0.08,'SubspaceDimension',2*nmode_fluid+25,'MaxIterations',400);
[freq_fluid,renumf] = sort(sqrt(diag(real(valp_f)))/(2*pi));
freq_fluid=real(freq_fluid); %natural acoustic frequencies
Phi_f = real(vecp_f); %eigenvectors
disp([freq_fluid(1:9)])
%----PLOT FLUID MODES
for i = 1:20 %number of fluid modes to plot
    matlab2VTK_fluid(fmesh.nodes,fmesh.elems,Phi_f(:,i),['Acoustic_Mode/Modefluid_',num2str(i)])
end
%%
%---Generate INTERFACE
%====================================================
%------------FSi surface from fluid mesh------------%
%====================================================
[nodesftheta,nodesfr]=cart2pol(nodesf(:,1),nodesf(:,2));
numfsif = find( (nodesfr(:) > Rint-1e-7 &  nodesfr(:) < Rint+1e-7));
logical_index_f = all(ismember(elems_skinf, numfsif), 2);
elemsfsif = elems_skinf(logical_index_f, :);
%======================================================
% Mesh fsi (from fluid)
%======================================================
nnodesf = size(nodesf,1);
ndoff=  size(nodesf,1);
Uf = zeros(3*nnodesf,1);
Pf = zeros(ndoff,1);
%---PLOT INTERFACE 1 (FROM FLUID)
matlab2VTK_interface(nodesf,elemsfsif,Uf,Pf,'MESH/fsimeshfluid') 
%======================================================
%------------FSi surface from structure mesh
%======================================================
[nodesstheta,nodessr]=cart2pol(nodess(:,1),nodess(:,2));
numfsis = find( (nodessr(:) > Rint-1e-7 &  nodessr(:) < Rint+1e-7));
logical_index_s = all(ismember(elems_skins, numfsis), 2);
elemsfsis = elems_skins(logical_index_s, :);
%======================================================
% Mesh fsi (from structure)
%======================================================
ndofs=  3*size(nodess,1);
Us = zeros(ndofs,1);
nnodess = size(nodess,1);
Ps = zeros(nnodess,1);
%---PLOT INTERFACE 2 (FROM STRUCTURE)
matlab2VTK_interface(nodess,elemsfsis,Us,Ps,'MESH/fsimeshstructure')
%--------------BUILD COUPLING MATRICES
[Cfsi,Area]=couplingMatrix3DMatching(elemsfsif,elemsfsis,ndoff,ndofs,nodess,nodesf,rhof);
%%
% VIBROACOUSTIC PROBLEM
%==============================%
%         (U,P) Formulation    %
Kz = sparse(zeros(size(K)));
dimensione_H = size(H);
Hz = spalloc(dimensione_H(1), dimensione_H(2), 0);
Cz = sparse(zeros(ndofs,ndoff));
dofp = [1:ndofs];
%---VIBROACOUSTIC PROBLEM
Mh=  [ M(dofu,dofu)  , Cz(dofu,:) ; 
      -Cfsi(dofu,:)' , Q(:,:)];
Kh = [K(dofu,dofu)   , Cfsi(dofu,:) ; 
      Cz(dofu,:)'    , H(:,:)];
%-----Normalization of structural modes
Vecps_normalized = zeros(ndofs,ndofs);
for i = 1:nmode_solid
    Vnorm1 = sqrt(real(Phi_s(:,i))'*M*real(Phi_s(:,i)));
    Vecps_normalized(:,i) = (1/Vnorm1)*real(Phi_s(:,i));
end
%-----Normalization of cavity modes
Vecpf_normalized = zeros(ndoff,ndoff);
for i = 1:nmode_fluid
    Vnormf = sqrt(real(Phi_f(:,i))'*Q*real(Phi_f(:,i)));
    Vecpf_normalized(:,i) = (1/Vnormf)*real(Phi_f(:,i));
end
%%
%-------------------------------------------%
%-----------REDUCED ORDER MODEL-------------%
%-------------------------------------------%
S = 21; %STRUCTURAL MODES TAKE INTO ACCOUNT FOR REDUCED PROBLEM
F = 115;%ACOUSTIC MODES TAKE INTO ACCOUNT FOR REDUCED PROBLEM
Phi_s_c = Vecps_normalized(:,1:S); 
Phi_f_c = Vecpf_normalized(:,1:F);
Kr    = diag(diag(Phi_s_c'*K*Phi_s_c)); % REDUCED STRU. STIFF. MAT.
Mr    = diag(diag(Phi_s_c'*M*Phi_s_c)); % REDUCED STRU. MASS MAT
Hr    = diag(diag(Phi_f_c'*H*Phi_f_c)); % REDUCED FLUID STIFF. MAT
Qr    = diag(diag(Phi_f_c'*Q*Phi_f_c)); % REDUCED FLUID MASS MAT
CfsrK = Phi_s_c'*Cfsi*Phi_f_c;          % REDUCED COUPLING MAT.
CfsrM = Phi_f_c'*Cfsi'*Phi_s_c;         % REDUCED COUPLING MAT.
CzrK  = Phi_f_c'*Cz'*Phi_s_c;
CzrM  = Phi_s_c'*Cz*Phi_f_c;
%------REDUCED VIBROACOUSTIC MATRICES-------%
MhROM= [Mr      , CzrM ; 
       -CfsrM    ,Qr];
KhROM=[Kr       , CfsrK ; 
       CzrK     , Hr];
%------REDUCED MODAL ANALYSIS-----%
[Vecpcoupl_red,valpcoupl_red] = eig(full(KhROM),full(MhROM));
[freq_coupl_red, renumc_red]=sort(sqrt(diag(valpcoupl_red))/(2*pi));
freq_coupl_red = real(freq_coupl_red);
disp([freq_coupl_red(1:20)]) %coupled modes
%---------------------------------------------------
% Modal analysis
Umatc_red = zeros(S,S);
Pmatc_red = zeros(F,F);
for i=1:length(freq_coupl_red)
    j = renumc_red(i);
    Umatc_red(:,i) = real(Vecpcoupl_red(1:S,j));
    Pmatc_red(:,i) = real(Vecpcoupl_red(S+1:S+F,j));
end
Umatc_red = Phi_s_c*Umatc_red;
Pmatc_red = Phi_f_c*Pmatc_red;
% =====================================================%
%          Generate paraview files for coupled modes   %
% =====================================================%
for i=1:30
    matlab2VTK_fluid(nodesf,elemsf,Pmatc_red(:,i),['Coupled_Mode_Reduced/ModeFluid',num2str(i)])
    matlab2VTK_structure(nodess,elemss,Umatc_red(:,i),['Coupled_Mode_Reduced/ModeSolid',num2str(i)])
end