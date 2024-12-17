%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-----------Isogeometric analysis Vibroacoustic problem-----------%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc
clear variables
close all
%%
addpath ../C_files/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../function_TL_CH/
%% 
%--------------STRUCTURE--------------%
%----Geometry Charateristics
R_s= 1;    %external radius
L  = 3;    %lenght
t  = 0.05; %thickness
%%
%-----CONSTRUCT GEOMETRY STRUCTURE AND PLOT INITIAL GEOMETRY
[NURBS_struct_INIT] = generateNURBSDataHollowCylinder(R_s,t,L);
generateIGAvisu(NURBS_struct_INIT,'Geometry/Structure',1);
%%
%---P REFINEMENT (ORIGINAL ORDER p=2,r=2 and r=2) 
%THE FINAL POLYNOMIAL ORDER WILL BE THE (ORIGINAL ORDER+ELE_ORD_)
clear NURBSnew_struct
ele_ord_p = 3; %polynomial order in the RADIAL DIRECTION 
ele_ord_q = 1; %polynomial order in the THICKNESS 
ele_ord_r = 3; %polynomial order in the AXIAL DIRECTION 
[NURBS_struct] = pRefine3D(NURBS_struct_INIT,ele_ord_p,ele_ord_q,ele_ord_r);

% 3D H-REFINEMENT
noGPs =  max([ele_ord_p,ele_ord_q,ele_ord_r])+2+1; %p+1
nex =  2 ;  % NUMBER OF ELEMENTS (radial per each volume)
ney =  2 ;  % NUMBER OF ELEMENTS (thickness)
nez =  2 ;  % NUMBER OF ELEMENTS (axial)
ref_ax   =   nez-1; %if you want 4 elemnt u need to add (4-1) knots
ref_circ =   nex-1;
ref_rad  =   ney-1;
for i = 1:length(NURBS_struct) %REFINEMENT FUNCTION
    [NURBSnew_struct(i)]=refine3D(NURBS_struct(i),ref_circ,ref_rad,ref_ax);
end
%----PLOT Refined Geometry
generateIGAvisu(NURBSnew_struct,'Mesh\Refinement\Structure',1);
%%
%--------STRUCTURAL PROPERTIES-------%
E   = 7e10; %YOUNG MODULUS
nu  = 0.3;  %POISSON COEFFICIENT
rho = 2700; %STRUCTURAL DENSITY
lamb = nu*E/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));
C_str =[lamb+2*mu lamb      lamb      0  0  0;
        lamb      lamb+2*mu lamb      0  0  0;
        lamb      lamb      lamb+2*mu 0  0  0;
        0         0         0         mu 0  0;
        0         0         0         0  mu 0;
        0         0         0         0  0  mu];
%-Build Elementary Structural Matrices
for i = 1:length(NURBSnew_struct)
    % KM_matrices_3D
    [noDofs,noElems,Ks,Ms] = KM_matrices_3D_parallel(NURBSnew_struct(i),noGPs,C_str,rho);    
    NURBSnew_struct(i).K = Ks;
    NURBSnew_struct(i).M = Ms;
    NURBSnew_struct(i).noCtrPts = 3*size(NURBSnew_struct(i).controlPts,1);
end
%------ASSEMBLE GLOBAL STRUCTURAL MATRICES------%
[IGA_s,NURBSnew_struct] = AssembleStructMatrixIGA(NURBSnew_struct);
doftot_s=[1:IGA_s.nbB];
%%
%------BOUNDARY CONDITIONS (x,y)
%--nodes definition order(x1,x2,x3,...xn,y1,y2,y3,yn,...z1,z2,z3,z4,...zn)
z0= 0;
zL= L;
tol=10^-6;
dofb0 = find(ismembertol(IGA_s.B(:,3),z0,tol));
dofbL = find(ismembertol(IGA_s.B(:,3),zL,tol));
dofbx = vertcat(dofb0,dofbL);
dofby = vertcat(dofb0+size(IGA_s.B,1),dofbL+size(IGA_s.B,1));
dofb = vertcat(dofbx,dofby); %blocked dofs
dofu = setdiff(doftot_s,dofb); %free dofs
%%
%===============================================%
% Compute eigenvalue structure analysis
%===============================================%
nbmode_solid = 10;
[Vecp_solid,Valp_solid] = eigs(IGA_s.K(dofu,dofu),IGA_s.M(dofu,dofu),nbmode_solid,'smallestreal','tolerance',0.08,'SubspaceDimension',2*nbmode_solid+60,'MaxIterations',600);
[freq_solid,renums] = sort(sqrt(diag(real(Valp_solid)))/(2*pi));
freq_solid=real(freq_solid); %natural frequencies
disp([freq_solid(1:10)]);
%===============================================%
%                Structural modes               %
%===============================================%
PhiStructMat = zeros(IGA_s.nbB,nbmode_solid);
PhiStructMat(dofu,:) = real(Vecp_solid);
%---------------------------------------------------
Umats = zeros(IGA_s.nbB,nbmode_solid);
for i=1:nbmode_solid
    j = renums(i);
    Umats(:,i) = PhiStructMat(:,j);
end
%-----PLOT STRUCTURAL MODES
for j=1:10  %number of modes to genrate in paraview
    for i = 1:length(NURBSnew_struct)
        ndof = 3*size(NURBSnew_struct(i).controlPts,1);
        NURBSnew_struct(i).ncontrolPts = size(NURBSnew_struct(i).controlPts,1);
        Ux = real(Umats(NURBSnew_struct(i).global_to_local,j));
        Uy = real(Umats(NURBSnew_struct(i).global_to_local+size(IGA_s.B,1),j));
        Uz = real(Umats(NURBSnew_struct(i).global_to_local+2*size(IGA_s.B,1),j));
        U  = reshape([Ux,Uy,Uz]',[ndof,1]);
        NURBSnew_struct(i).Ux = Ux;
        NURBSnew_struct(i).Uy = Uy;
        NURBSnew_struct(i).Uz = Uz;
        NURBSnew_struct(i).U  = U;
    end
    generateIGAvisuField(NURBSnew_struct,'StructuralMode\Refinement',j);
end
%%
%--------------CAVITY--------------%
%----Geometry Charateristics
R_f    = R_s-t; %Radius
L      = 3;
clear NURBSnew_fluid
%-----CONSTRUCT GEOMETRY CAVITY  AND PLOT GEOMETRY
[NURBS_fluid_INIT] = generateNURBSDataCylinder(R_f,L);
generateIGAvisu(NURBS_fluid_INIT,'Geometry/Cavity',1); %Visualize cavity
%---P REFINEMENT (ORIGINAL ORDER p=2,r=2 and r=2)
ele_ord_p_f = ele_ord_p; %RADIAL DIRECTION
ele_ord_q_f = ele_ord_p; %THICKNESS
ele_ord_r_f = ele_ord_p; %AXIAL
[NURBS_fluid] = pRefine3D(NURBS_fluid_INIT,ele_ord_p_f,ele_ord_q_f,ele_ord_r_f);
% 3D H-REFINEMENT
for i = 1:length(NURBS_fluid)
   [NURBSnew_fluid(i),Tnew]=refine3D(NURBS_fluid(i),ref_circ,ref_circ,ref_ax);
end
%--PLOT REFINED Geometry
generateIGAvisu(NURBSnew_fluid,'Mesh/Refinement/Cavity',2);
%%
%==========BUILD ACOUSTIC MATRICES===========%%%
c0   = 340; %SPEED SOUND (m/s)
rhof = 1.2; %fluid density (kg/m^3)
%------Build the elementary Acoustic Matrices
for i = 1:length(NURBSnew_fluid)
        [~,Ka,Ma] = HQ_matrices_3D_parallel(NURBSnew_fluid(i),noGPs,c0);
       NURBSnew_fluid(i).K        = Ka;
       NURBSnew_fluid(i).M        = Ma;
       NURBSnew_fluid(i).noCtrPts = size(NURBSnew_fluid(i).controlPts,1);
end
%------ASSEMBLE GLOBAL ACOUSTIC MATRICES------%
[IGA_f,NURBSnew_fluid] = AssembleFluidMatrixIGA(NURBSnew_fluid);
% ===============================================%
% Compute eigenvalue cavity analysis
% ===============================================%
nbmode_fluid = 10;
[vecp_f,valp_f] = eigs(IGA_f.K,IGA_f.M,nbmode_fluid,'smallestreal','tolerance',0.08,'MaxIterations',600);
[freq_fluid,renumf] = sort(sqrt(diag(real(valp_f)))/(2*pi));
freq_fluid=real(freq_fluid);
disp([freq_fluid(1:9)])
%==================================================%
%                  Visualize fluid MODES           %
%==================================================%
PhifFluidMat = zeros(IGA_f.nbB,nbmode_fluid);
PhifFluidMat(:,:) = real(vecp_f(:,1:nbmode_fluid));
%---------------------------------------------------
Umatf = zeros(IGA_f.nbB,nbmode_fluid);
for i=1:nbmode_fluid
    j = renumf(i);
    Umatf(:,i) = PhifFluidMat(1:IGA_f.nbB,j);
end
Umatf=real(Umatf);
% ----PLOT FLUID MODES
for j=1:10
    for i = 1:length(NURBSnew_fluid)
        ndof = 3*size(NURBSnew_fluid(i).controlPts,1);
        NURBSnew_fluid(i).ncontrolPts = size(NURBSnew_fluid(i).controlPts,1);
        Px = real(Umatf(NURBSnew_fluid(i).global_to_local,j));
        Py = zeros(length(Px),1);
        Pz = zeros(length(Px),1);
        P  = reshape([Px,Py,Pz]',[ndof,1]);
        NURBSnew_fluid(i).Ux = Px;
        NURBSnew_fluid(i).Uy = Py;
        NURBSnew_fluid(i).Uz = Pz;
        NURBSnew_fluid(i).U  = P;
    end
    generateIGAvisuField(NURBSnew_fluid,'AcousticMode/Refinement',j);
end
%%
%------COUPLING INTERFACE------%
% EXTRACTION SURFACES FROM 3D GEOMETRY
[NURBSnew_fluid] = extractNURBS2D(NURBSnew_fluid);
% Selection the surfaces involved in the coupling
% the first coloumn is the volume involed and
% the second coloumn is the surface involved
selection  = [1, 4;
              2, 4;
              3, 4;
              4, 4];
INTERFACE = [];
for k = 1:size(selection,1)
    i = selection(k,1);
    j = selection(k,2);
    INTERFACE(k).NURBS2D = NURBSnew_fluid(i).NURBS2D(j);
    INTERFACE(k).volinfo = i;
    INTERFACE(k).surfinfo= j;
end
%---PLOT INTERFACE
generateIGAvisu2D(INTERFACE,'Mesh/Refinement/Interface',1,noGPs)
%%
%-----C MATRIX FOR EACH VOLUME
%---BUILD ELEMENTARY COUPLING MATRICES
for i=1:size(INTERFACE,2)
    [Cfs,PhyGP,Normals,Area,idF,idS] = C_matrix_curved(INTERFACE(i).NURBS2D.controlPts, ...
                                                       size(INTERFACE(i).NURBS2D.controlPts,1), ...
                                                       size(INTERFACE(i).NURBS2D.controlPts,1), ...
                                                       INTERFACE(i).NURBS2D.controlPts, ...
                                                       INTERFACE(i).NURBS2D.controlPts, ...
                                                       INTERFACE(i).NURBS2D.uknot, ...
                                                       INTERFACE(i).NURBS2D.vknot, ...
                                                       INTERFACE(i).NURBS2D.noPtsX, ...
                                                       INTERFACE(i).NURBS2D.noPtsY, ...
                                                       INTERFACE(i).NURBS2D.p,INTERFACE(i).NURBS2D.q,INTERFACE(i).NURBS2D.weights,noGPs);
    INTERFACE(i).NURBS2D.Cfs = Cfs;
    INTERFACE(i).NURBS2D.idF = idF;
    INTERFACE(i).NURBS2D.Area = Area;
    INTERFACE(i).NURBS2D.PhyGP = PhyGP;
    INTERFACE(i).NURBS2D.Normals = Normals;
end
%----ASSAMBLE COUPLING MATRIX
for i = 1:size(INTERFACE,2)
    [IGA_c] = AssembleCouplingIGA(INTERFACE(i),IGA_s,IGA_f);
end
% ==========================================================%
%             ASSEMBLE VIBROACOUSTIC PROBLEM                %
% ==========================================================%
Kz = sparse(zeros(size(IGA_s.K)));
dimensione_H = size(IGA_f.K);
Hz = spalloc(dimensione_H(1), dimensione_H(2), 0);
Cz = sparse(zeros(IGA_s.nbB,IGA_f.nbB));
Mh=  [ IGA_s.M(dofu,dofu)     , Cz(dofu,:) ; 
       -IGA_c.Cfs(dofu,:)'    , (1/rhof)*IGA_f.M(:,:)];
   
Kh = [ IGA_s.K(dofu,dofu)     , IGA_c.Cfs(dofu,:) ; 
       Cz(dofu,:)'            , (1/rhof)*IGA_f.K(:,:)];
%-----Normalization of structural modes
Vecps_normalized = zeros(IGA_s.nbB,IGA_s.nbB);
for i = 1:nbmode_solid
    Vnorm1 = sqrt(real(PhiStructMat(:,i))'*IGA_s.M*real(PhiStructMat(:,i)));
    Vecps_normalized(:,i) = (1/Vnorm1)*real(PhiStructMat(:,i));
end
%-----Normalization of cavity modes
Vecpf_normalized = zeros(IGA_f.nbB,IGA_f.nbB);
for i = 1:nbmode_fluid
    Vnormf = sqrt(real(PhifFluidMat(:,i))'*IGA_f.M*real(PhifFluidMat(:,i)));
    Vecpf_normalized(:,i) = (1/Vnormf)*real(PhifFluidMat(:,i));
end
%%
%%-------------REDUCED ORDER MODELS---------------%%
%---------CONSTRUCTION OF REDUCED MATRICES---------%
S = 50;  %STRUCTURAL MODES TO TAKE INTO ACCOUNT IN THE REDUCED PROBLEM
F = 50; %FLUID MODES TO TAKE INTO ACCOUNT IN THE REDUCED PROBLEM
Phi_s_c = Vecps_normalized(:,1:S); 
Phi_f_c = Vecpf_normalized(:,1:F);
Kr    = diag(diag(Phi_s_c'*IGA_s.K*Phi_s_c)); %REDUCED STRUCT. STIFF. MAT.
Mr    = diag(diag(Phi_s_c'*IGA_s.M*Phi_s_c)); %REDUCED STRUCT. MASS MAT.
Hr    = diag(diag(Phi_f_c'*IGA_f.K*Phi_f_c)); %REDUCED ACOU.   STIFF. MAT.
Qr    = diag(diag(Phi_f_c'*IGA_f.M*Phi_f_c)); %REDUCED ACOU.   MASS MAT.
CfsrK = Phi_s_c'*IGA_c.Cfs*Phi_f_c;           %REDUCED COUPL. MAT.
CfsrM = Phi_f_c'*IGA_c.Cfs'*Phi_s_c;          %REDUCED COUPL. MAT.
CzrK  = sparse(Phi_f_c'*Cz'*Phi_s_c);         %REDUCED COUPL. MAT.
CzrM  = sparse(Phi_s_c'*Cz*Phi_f_c);          %REDUCED COUPL. MAT.
Kzr   = sparse(zeros(size(Kr)));
Hzr   = sparse(zeros(size(Hr)));
%--REDUCED SYSTE,
MhROM= [Mr      ,          CzrM ; 
        -CfsrM   , (1/rhof)*Qr];
KhROM=[Kr       ,        CfsrK ; 
       CzrK     , (1/rhof)*Hr];
MhROM=sparse(MhROM);
KhROM=sparse(KhROM);
%----MODAL ANALYSIS OF REDUCED PROBLEM
[Vecpcoupl_red,valpcoupl_red] = eig(full(KhROM),full(MhROM));
[freq_coupl_red, renumc_red]=sort(sqrt(diag(valpcoupl_red))/(2*pi));
freq_coupl_red = real(freq_coupl_red);
disp([freq_coupl_red(1:20)]);
%PLOT REDUCED MODES
%---------------------------------------------------
% Modal analysis
Pmatc_red = zeros(F,F);
Umatc_red = zeros(S,S);
for i=1:length(freq_coupl_red)
    j = renumc_red(i);
    Umatc_red(:,i) = real(Vecpcoupl_red(1:S,j));
    Pmatc_red(:,i) = real(Vecpcoupl_red(S+1:S+F,j));
end
Umatc_red = Phi_s_c*Umatc_red;
Pmatc_red = Phi_f_c*Pmatc_red;
%%
%=======================================================%
%               Visualize COUPLED MODES                 %
%=======================================================%
%----------STRUCTURAL SIDE
for j = 1:10
    for i = 1:length(NURBSnew_struct)
            ndof = 3*size(NURBSnew_struct(i).controlPts,1);
            NURBSnew_struct(i).ncontrolPts = size(NURBSnew_struct(i).controlPts,1);
            Uxc = real(Umatc_red(NURBSnew_struct(i).global_to_local,j));
            Uyc = real(Umatc_red(NURBSnew_struct(i).global_to_local+size(IGA_s.B,1),j));
            Uzc = real(Umatc_red(NURBSnew_struct(i).global_to_local+2*size(IGA_s.B,1),j));
            Uc  = reshape([Uxc,Uyc,Uzc]',[ndof,1]);
            NURBSnew_struct(i).Ux = Uxc;
            NURBSnew_struct(i).Uy = Uyc;
            NURBSnew_struct(i).Uz = Uzc;
            NURBSnew_struct(i).U  = Uc;
    end
    generateIGAvisuField(NURBSnew_struct,'Coupled_Mode_Reduced/Refinement/ModeStructure',j);
end
% ----------FLUID SIDE
for j = 1:10
    for i = 1:length(NURBSnew_fluid)
        ndof = 3*size(NURBSnew_fluid(i).controlPts,1);
        NURBSnew_fluid(i).ncontrolPts = size(NURBSnew_fluid(i).controlPts,1);
        Pxc = real(Pmatc_red(NURBSnew_fluid(i).global_to_local,j));
        Pyc = zeros(length(Pxc),1);
        Pzc = zeros(length(Pxc),1);
        Pc  = reshape([Pxc,Pyc,Pzc]',[ndof,1]);
        NURBSnew_fluid(i).Ux = Pxc;
        NURBSnew_fluid(i).Uy = Pyc;
        NURBSnew_fluid(i).Uz = Pzc;
        NURBSnew_fluid(i).U  = Pc;
    end
    generateIGAvisuField(NURBSnew_fluid,'Coupled_Mode_Reduced/Refinement/ModeFluid',j);
end