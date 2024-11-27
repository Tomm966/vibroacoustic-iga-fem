function [nodes, elems] = generateQuadMesh(Lx, Ly, nelemx, nelemy)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% Lx: Length of the mesh in the x-direction.
% Ly: Length of the mesh in the y-direction.
% nelemx: Number of elements along the x-direction.
% nelemy: Number of elements along the y-direction.
% OUTPUT:
% nodes: An N×2 matrix containing the coordinates of the nodes, where N is the total number of nodes in the mesh.
% elems: An M×4 matrix containing the connectivity of the elements, where M is the total number of quadrilateral elements.

    % Generate a 2D mesh with quad 4 elements.

    % Create nodes
    dx = Lx / nelemx;
    dy = Ly / nelemy;

    x = 0:dx:Lx;
    y = 0:dy:Ly;

    [X, Y] = meshgrid(x, y);
    nodes = [X(:), Y(:)];

    % Create elements
    elems = zeros(nelemx * nelemy, 4);

    for j = 1:nelemx
        for i = 1:nelemy
            n1 = (j - 1) * (nelemy + 1) + i;
            n2 = n1 + 1;
            n3 = n1 + nelemy + 1;
            n4 = n3 + 1;

            elems((i - 1) * nelemx + j, :) = [n1, n3, n4, n2];
        end
    end
end