function [] = generateIGAvisuField(NURBS,namefolder,num)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% NURBS: An array of structures, each containing data for a NURBS object, including control points, number of points in each direction (X, Y, Z), and other parameters related to the geometry and field.
% namefolder: The name of the output directory where VTK files will be saved.
% num: An identifier for naming the output files.
% OUTPUT:
% The function generates VTK files for:
% Control net visualization.
% Patch edges.
% Patch volumes.
% Additional mesh data (commented out).


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
    UAll          = [];
    for i = 1:length(NURBS)
        elemsnet = NURBS(i).elemsnet; %check if it crash somehow
        elemsnetAll   = [elemsnetAll ; NURBS(i).elemsnet + size(controlPtsAll,1)];
        controlPtsAll = [controlPtsAll ; NURBS(i).controlPts];
        PIDAll        = [PIDAll; i*ones(size(NURBS(i).elemsnet,1),1)];
        UAll          = [UAll  ;NURBS(i).U];  
    end
    matlab2VTK_net(controlPtsAll,elemsnetAll,UAll,PIDAll,['./',namefolder,'/net_patch',num2str(num)])
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
    [nodes_param,elems_param,elems_skin_param] = generateFEMeshQuadratique20(1.0,1.0,1.0,15,15,15);
    nodesvAll = [];
    fieldvAll = [];
    elemsvAll = [];
    PIDvAll   = [];
    for i = 1:length(NURBS)
        [nodes_vol] = createPatchVolume(NURBS(i),nodes_param);
        [field_vol] = createPatchField(NURBS(i),nodes_param,NURBS(i).Ux,NURBS(i).Uy,NURBS(i).Uz);
        elemsvAll = [elemsvAll ; elems_param + size(nodesvAll,1)];
        nodesvAll = [nodesvAll ; nodes_vol];
        fieldvAll = [fieldvAll ; field_vol];
        PIDvAll   = [PIDvAll; i*ones(size(elems_param,1),1)];
    end
    U = reshape([fieldvAll(:,1),fieldvAll(:,2),fieldvAll(:,3)]',[3*size(fieldvAll,1),1]);
    matlab2VTK_Carrera(nodesvAll(:,[1,2,3]),elemsvAll,...
                         U,...
                         zeros(3*size(nodesvAll,1),1),...
                         PIDvAll,['./',namefolder,'/nurbs',num2str(num)]);   
    %======================================
    % Generate elems edges
    %======================================
%     nodeseAll = [];
%     fieledAll = [];
%     elemseAll = [];
%     PIDeAll   = [];
%     for i = 1:length(NURBS)
%         [nodes_physic_e,elems_e]=generateElemsIGA(NURBS(i));
%         [field_physic_e,elems_e]=generateElemsIGAField(NURBS(i),NURBS(i).Ux,NURBS(i).Uy,NURBS(i).Uz);
%         elemseAll = [elemseAll ; elems_e + size(nodeseAll,1)];
%         nodeseAll = [nodeseAll ; nodes_physic_e];
%         fieledAll = [fieledAll ; field_physic_e];
%         PIDeAll   = [PIDeAll; i*ones(size(elems_e,1),1)];
%     end
%     U = reshape([fieledAll(:,1),fieledAll(:,2),fieledAll(:,3)]',[3*size(fieledAll,1),1]);
%     matlab2VTK_net(nodeseAll(:,[1,2,3]),...
%                    elemseAll,...
%                    U,...
%                    ones(size(elemseAll,1),1),['./',namefolder,'/elems',num2str(num)])
end