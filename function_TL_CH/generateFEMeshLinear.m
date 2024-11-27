function [nodes,elems,elems_skin] = generateFEMeshLinear(Lx,Ly,Lz,nex,ney,nez)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% Lx: Length of the mesh in the x-direction.
% Ly: Length of the mesh in the y-direction.
% Lz: Length of the mesh in the z-direction.
% nex: Number of elements along the x-axis.
% ney: Number of elements along the y-axis.
% nez: Number of elements along the z-axis.

% OUTPUT:
% nodes: Matrix of node coordinates (size: nnodes x 3), where nnodes = (nex+1) * (ney+1) * (nez+1) and each row represents a node in 3D space (x, y, z).
% elems: Element connectivity matrix (size: nelem x 8), where nelem = nex * ney * nez and each row contains the 8 node indices that define each hexahedral element.
% elems_skin: Matrix of surface (skin) elements (size: n_skin_elements x 4), where each row contains the 4 node indices that define the quadrilateral surface elements for the boundaries of the mesh.
%-------------------------------
    nnodex = nex+1;
    nnodey = ney+1;
    nnodez = nez+1;
    %-------------------------------
    Xi     = linspace(0,Lx,nnodex);
    Yi     = linspace(0,Ly,nnodey);
    Zi     = linspace(0,Lz,nnodez);
    %-------------------------------
    nnodes = nnodex*nnodey*nnodez;
    %===============================
    % Nodes coordinates
    %===============================
    A = 1;
    nodes  = zeros(nnodes,3);
    for iz = 1:nnodez
        for iy = 1:nnodey
            for ix = 1:nnodex
                nodes(A,1) = Xi(ix);
                nodes(A,2) = Yi(iy);
                nodes(A,3) = Zi(iz);
                A = A+1;
            end
        end
    end
    %===============================
    % Element Connectivity
    %===============================
    nelem = nex*ney*nez;
    elems = zeros(nelem,8);
    elem_count = 1;
    for iz = 1:nez
        for iy = 1:ney
            for ix = 1:nex
                n1 = (iz-1)*nnodey*nnodex + (iy-1)*nnodex + ix;
                n2 = n1 + 1;
                n3 = n2 + nnodex;
                n4 = n1 + nnodex;
                n5 = n1 + nnodey*nnodex;
                n6 = n2 + nnodey*nnodex;
                n7 = n3 + nnodey*nnodex;
                n8 = n4 + nnodey*nnodex;
                elems(elem_count,:) = [n1 n2 n3 n4 n5 n6 n7 n8];
                elem_count = elem_count + 1;
            end
        end
    end
    %===============================
    % Extract skin mesh
    %===============================
    %skin bottom
    nelem_bottom = nex*ney;
    elems_bottom = zeros(nelem_bottom,4);
    elem_count = 1;
    for iy = 1:ney
        for ix = 1:nex
            n1 = (iy-1)*nnodex + ix;
            n2 = n1+1;
            n3 = n2+nnodex;
            n4 = n1+nnodex;
            elems_bottom(elem_count,:) = [n1 n4 n3 n2];
            elem_count = elem_count + 1;
        end
    end
    %skin up
    nelem_up = nex*ney;
    elems_up = zeros(nelem_up,4);
    elem_count = 1;
    for iy = 1:ney
        for ix = 1:nex
            n1 = nez*nnodey*nnodex + (iy-1)*nnodex + ix;
            n2 = n1+1;
            n3 = n2+nnodex;
            n4 = n1+nnodex;
            elems_up(elem_count,:) = [n1 n2 n3 n4];
            elem_count = elem_count + 1;
        end
    end
    %skin left
    nelem_left = ney*nez;
    elems_left = zeros(nelem_left,4);
    elem_count = 1;
    for iz = 1:nez
        for iy = 1:ney
            n1 = (iy-1)*nnodex + 1 + (iz-1)*nnodey*nnodex;
            n2 = n1+ nnodex;
            n3 = n2+ nnodex*nnodey;
            n4 = n1+ nnodex*nnodey;
            elems_left(elem_count,:) = [n1 n4 n3 n2];
            elem_count = elem_count + 1;
        end
    end
    %skin front
    nelem_front = nex*nez;
    elems_front = zeros(nelem_front,4);
    elem_count = 1;
    for iz = 1:nez
        for ix = 1:nex
            n1 = ix + (iz-1)*nnodex*nnodey;
            n2 = n1 + 1;
            n3 = n2 + nnodex*nnodey;
            n4 = n1 + nnodex*nnodey;
            elems_front(elem_count,:) = [n1 n2 n3 n4];
            elem_count = elem_count + 1;
        end
    end
    %skin right
    nelem_right = ney*nez;
    elems_right = zeros(nelem_right,4);
    elem_count = 1;
    for iz = 1:nez
        for iy = 1:ney
            n1 = iy*nnodex + (iz-1)*nnodey*nnodex;
            n2 = n1+ nnodex;
            n3 = n2+ nnodex*nnodey;
            n4 = n1+ nnodex*nnodey;
            elems_right(elem_count,:) = [n1 n2 n3 n4];
            elem_count = elem_count + 1;
        end
    end
    %skin back
    nelem_back = nex*nez;
    elems_back = zeros(nelem_back,4);
    elem_count = 1;
    for iz = 1:nez
        for ix = 1:nex
            n1 = ix + iz*nnodex*nnodey - nnodex;
            n2 = n1 + 1;
            n3 = n2 + nnodex*nnodey;
            n4 = n1 + nnodex*nnodey;
            elems_back(elem_count,:) = [n2 n1 n4 n3];
            elem_count = elem_count + 1;
        end
    end
    %total skin
    elems_skin = [elems_bottom;elems_up;elems_left;elems_front;elems_right;elems_back];
    
end