function [] = generateIGAvisu2D(INTERFACE,namefolder,num,noGPs)

% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% INTERFACE: An array of structures containing NURBS 2D data (control points, knot vectors, weights, etc.).
% namefolder: The name of the output directory where VTK files will be saved.
% num: An identifier for naming the output files.
% noGPs: The number of Gauss points for numerical integration.
% OUTPUT:
% The function does not return any outputs directly but generates VTK files for:
% Control net visualization.
% Surface meshes.
% Normals at Gauss points.
% Element edges.

    mkdir(namefolder)
    %======================================
    % Generate control net vtk all
    %======================================
    controlPtsAll = [];
    elemsnetAll   = [];
    PIDAll        = [];
    for j = 1:length(INTERFACE)
        elemsnetAll   = [elemsnetAll ; INTERFACE(j).NURBS2D.elemsnet + size(controlPtsAll,1)];
        controlPtsAll = [controlPtsAll ; INTERFACE(j).NURBS2D.controlPts];
        elemsnet      = INTERFACE(j).NURBS2D.elemsnet;
        PIDAll        = [PIDAll; j*ones(size(elemsnet,1),1)];
    
    end
    U = zeros(size(controlPtsAll,1)*3,1);
    matlab2VTK_net(controlPtsAll,elemsnetAll,U,PIDAll,['./',namefolder,'/net_patch',num2str(num)])
    %======================================
    % Generate surfaces all
    %======================================
    [nodes_param,elems_param] = generateQuadMesh(1.0,1.0,10,10);
    nodesvAll = [];
    elemsvAll = [];
    PIDvAll   = [];
    PvAll     = [];
    for j = 1:length(INTERFACE)
        [nodes_sur] = createPatchSurface(INTERFACE(j).NURBS2D,nodes_param);
        store = INTERFACE(j).NURBS2D.controlPts(:,1); %
        INTERFACE(j).NURBS2D.controlPts(:,1) = INTERFACE(j).NURBS2D.weights;
        [press_sur] = createPatchSurface(INTERFACE(j).NURBS2D,nodes_param);
        INTERFACE(j).NURBS2D.controlPts(:,1) = store ;%
        elemsvAll = [elemsvAll ; elems_param + size(nodesvAll,1)];
        nodesvAll = [nodesvAll ; nodes_sur];
        PIDvAll   = [PIDvAll; j*ones(size(elems_param,1),1)];
        PvAll   = [PvAll; press_sur(:,1)];
    end
    U = zeros(size(nodesvAll,1)*3,1);
    P = zeros(size(nodesvAll,1),1);
%     matlab2VTK_Quad4(nodesvAll(:,[1:3]),elemsvAll,U,PvAll, ...
%                 ['./',namefolder,'/nurbs',num2str(num)])%%
    matlab2VTK_Quad4(nodesvAll(:,[1:3]),elemsvAll,U,PvAll, ...
                ['./',namefolder,'/nurbs',num2str(num)])%%
    %======================================
    % External normal at Gauss Points
    %======================================
    area = 0;
    GPall = [];
    Uall  = [];
    for j = 1:length(INTERFACE)
        NURBS2D = INTERFACE(j).NURBS2D;
        [PhyGP, Normals,Area] = Area_curved(NURBS2D.controlPts, ...
                                            NURBS2D.uknot,  ...
                                            NURBS2D.vknot,  ...
                                            NURBS2D.noPtsX, ...
                                            NURBS2D.noPtsY, ...
                                            NURBS2D.p,NURBS2D.q, ...
                                            NURBS2D.weights,noGPs);
    
        NURBS2D.Area = Area;
        area  =area + Area;
        NURBS2D.PhyGP = PhyGP;
        NURBS2D.Normals = Normals;
        U = reshape(Normals,3*size(PhyGP',1),1);
        INTERFACE(j).NURBS2D = NURBS2D;
        GPall = [GPall;PhyGP'];
        Uall  = [Uall;U];
    end
    matlab2VTK_points(GPall,[1:size(GPall,1)]',Uall, ...
            ['./',namefolder,'/normal_GP',num2str(num)])
    %======================================
    % Generate elems edges
    %======================================
    nodeseAll = [];
    elemseAll = [];
    PIDeAll   = [];
    for i = 1:length(INTERFACE)
        NURBS2D = INTERFACE(i).NURBS2D;
        [nodes_physic_e,elems_e] = generateElemsIGA2d(NURBS2D.controlPts, ...
                                            NURBS2D.uknot,  ...
                                            NURBS2D.vknot,  ...
                                            NURBS2D.noPtsX, ...
                                            NURBS2D.noPtsY, ...
                                            NURBS2D.p,NURBS2D.q, ...
                                            NURBS2D.weights);
        %[field_physic_e,elems_e]=generateElemsIGAField(NURBS(i),NURBS(i).Ux,NURBS(i).Uy,NURBS(i).Uz);
        elemseAll = [elemseAll ; elems_e + size(nodeseAll,1)];
        nodeseAll = [nodeseAll ; nodes_physic_e];
        %fieledAll = [fieledAll ; field_physic_e];
        PIDeAll   = [PIDeAll; i*ones(size(elems_e,1),1)];
    end
    U = zeros(3*size(nodeseAll,1),1);
    matlab2VTK_net(nodeseAll(:,[1,2,3]),...
                   elemseAll,...
                   U,...
                   ones(size(elemseAll,1),1),['./',namefolder,'/elems',num2str(num)])
    
end