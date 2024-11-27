function [] = generateIGAvisu(NURBS,namefolder,num)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% NURBS: An array of NURBS structures, each containing properties such as control points and mesh parameters.
% namefolder: The name of the folder where the generated VTK files will be stored.
% num: A numeric identifier for the files (used for naming purposes).
% OUTPUT:
% The function does not return any outputs directly; instead, it generates and saves files to the specified folder.

    mkdir(namefolder)
    %======================================
    % Generate control net meshes
    %======================================
    for i = 1:length(NURBS)
        noPtsX = NURBS(i).noPtsX;
        noPtsY = NURBS(i).noPtsY;
        noPtsZ = NURBS(i).noPtsZ;
        [elemsnet] = createNetVisu(noPtsX,noPtsY,noPtsZ);
        NURBS(i).elemsnet = elemsnet;
    end
    %======================================
    % Generate all control nets
    %======================================
    controlPtsAll = [];
    elemsnetAll   = [];
    PIDAll        = [];
    for i = 1:length(NURBS)
        elemsnetAll   = [elemsnetAll ; NURBS(i).elemsnet + size(controlPtsAll,1)];
        controlPtsAll = [controlPtsAll ; NURBS(i).controlPts];
        PIDAll        = [PIDAll; i*ones(size(elemsnet,1),1)];
       
    end
    U = zeros(size(controlPtsAll,1)*3,1);
    matlab2VTK_net(controlPtsAll,elemsnetAll,U,PIDAll,['./',namefolder,'/net_patch',num2str(num)])
    %======================================
    % Generate all patch edged
    %======================================
    nodespAll = [];
    elemspAll = [];
    PIDpAll   = [];
    for i = 1:length(NURBS)
        [nodes_param_test,nodes_physic,elems] = createPatchEdges(NURBS(i),10);
        elemspAll = [elemspAll ; elems + size(nodespAll,1)];
        nodespAll = [nodespAll ; nodes_physic];
        PIDpAll   = [PIDpAll; i*ones(size(elems,1),1)];
    end
    U = zeros(size(nodespAll,1)*3,1);
    matlab2VTK_net(nodespAll(:,[1,2,3]),elemspAll,U,PIDpAll,['./',namefolder,'/patch_edges',num2str(num)])
    %======================================
    % Generate patch volulmes
    %======================================
    [nodes_param,elems_param,elems_skin_param] = generateFEMeshQuadratique20(1.0,1.0,1.0,10,10,10);
    nodesvAll = [];
    elemsvAll = [];
    PIDvAll   = [];
    for i = 1:length(NURBS)
        [nodes_vol] = createPatchVolume(NURBS(i),nodes_param);
        elemsvAll = [elemsvAll ; elems_param + size(nodesvAll,1)];
        nodesvAll = [nodesvAll ; nodes_vol];
        PIDvAll   = [PIDvAll; i*ones(size(elems_param,1),1)];
    end
    matlab2VTK_Carrera(nodesvAll(:,[1,2,3]),elemsvAll,...
                         zeros(3*size(nodesvAll,1),1),...
                         zeros(3*size(nodesvAll,1),1),...
                         PIDvAll,['./',namefolder,'/nurbs',num2str(num)]);   
    %======================================
    % Generate elems edges
    %======================================
    nodeseAll = [];
    elemseAll = [];
    PIDeAll   = [];
    for i = 1:length(NURBS)
        [nodes_physic_e,elems_e]=generateElemsIGA(NURBS(i));
        elemseAll = [elemseAll ; elems_e + size(nodeseAll,1)];
        nodeseAll = [nodeseAll ; nodes_physic_e];
        PIDeAll   = [PIDeAll; i*ones(size(elems_e,1),1)];
    end
    matlab2VTK_net(nodeseAll(:,[1,2,3]),...
                   elemseAll,...
                   zeros(3*size(nodeseAll,1),1),...
                   ones(size(elemseAll,1),1),['./',namefolder,'/elems',num2str(num)])
end