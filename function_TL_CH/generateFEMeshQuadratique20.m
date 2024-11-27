function [nodes_unique,elems_unique,elems_skin_unique] = generateFEMeshQuadratique20(Lx,Ly,Lz,nex,ney,nez)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% Lx: Length of the mesh in the x-direction.
% Ly: Length of the mesh in the y-direction.
% Lz: Length of the mesh in the z-direction.
% nex: Number of elements along the x-axis.
% ney: Number of elements along the y-axis.
% nez: Number of elements along the z-axis.

% OUTPUT:
% nodes_unique: Matrix of unique node coordinates (size: n_unique_nodes x 3), where each row represents the coordinates of a unique quadratic node in 3D space (x, y, z).
% elems_unique: Element connectivity matrix (size: nelem x 20), where each row contains the 20 unique node indices that define each quadratic hexahedral element.
% elems_skin_unique: Matrix of surface (skin) elements (size: n_skin_elements x 8), where each row contains the 8 unique node indices that define the quadratic surface elements for the boundaries of the mesh.

%===============================
    % Creat linear mesh first
    %===============================
    [nodes,elems,elems_skin] = generateFEMeshLinear(Lx,Ly,Lz,nex,ney,nez);
    %===============================
    % From linear to quadratic
    %===============================
    nelem = size(elems,1);
    nnodes = size(nodes,1);
    new_nodes = zeros(12*nelem,3);
    elems_quadratic = zeros(nelem,20);
    for e = 1:nelem
        nums = [1 2;  %9
                2 3;  %10
                3 4;  %11
                4 1;  %12
                5 6;  %17
                6 7;  %18
                7 8;  %19
                8 5;  %20
                1 5;  %13
                2 6;  %14
                3 7;  %15
                4 8]; %16
        new_nodes((e-1)*12+[1:12],:) = 1/2*(nodes(elems(e,nums(:,1)),:) + ...
                                             nodes(elems(e,nums(:,2)),:));
        elems_quadratic(e,:) = [elems(e,:), nnodes+(e-1)*12+[1 2 3 4 9 10 11 12 5 6 7 8]];
%         keyboard
    end
    %===============================
    % New connectivity table with 
    % unique nodes
    %===============================
    nodes_quadratic = [nodes;new_nodes];
    [nodes_unique,~,renum] = uniquetol(nodes_quadratic,'ByRows',1e-8);
    elems_unique = renum(elems_quadratic);
    elems_skin_quadratic   = renum(elems_skin);
    %===============================
    % From linear to quadratic skin
    %===============================
    nelem_skin = size(elems_skin,1);
    elems_skin_unique = zeros(nelem_skin,8);
    for e = 1:nelem_skin
        nums = [1 2;   %5
                2 3;   %6
                3 4;   %7
                4 1];  %8
         [~,num] = ismembertol(1/2*(nodes_unique(elems_skin_quadratic(e,nums(:,1)),:) + ...
                                    nodes_unique(elems_skin_quadratic(e,nums(:,2)),:)), ...
                                    nodes_unique,"ByRows",1e-8);
         elems_skin_unique(e,:) = [elems_skin_quadratic(e,:), num'];
    end
       
end