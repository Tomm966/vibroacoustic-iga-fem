function [ output_args ] = matlab2VTK_interface(nodes,elems,U,P,filename)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% Inputs
% nodes: An n×3 matrix containing the x,y,z coordinates of the nodes in the mesh. Each row corresponds to a node.
% elems: An m×N matrix where each row represents an element by listing the indices of its nodes.
% U: A vector or matrix representing the displacements of each node. It is reshaped into a matrix of size n×3, where each row corresponds to the displacement in the x,y,z directions for each node.
% P: A vector of pressure values associated with each node, used for scalar visualization in ParaView.
% filename: The desired name of the output VTK file (without the .vtk extension).
 
% Outputs
% The function generates a .vtk file with the specified filename containing the mesh, displacement, and pressure data.

nnode = size(nodes,1);
U = reshape(U,[3,nnode])';

%=========================================================================%
%                    VISUALIZATION WITH PARAVIEW                          %
%=========================================================================%

%-------------------------------------------------------------------------%
%         Matlab read data from 1 to N                                    %
%         Paraview read data from 0 to N-1                                %
%-------------------------------------------------------------------------%

elems_para = elems - ones(size(elems)); % Renumerotation paraview

%-------------------------------------------------------------------------%
%                               Paraview quad8                            %
%-------------------------------------------------------------------------%
%                                                                         %
%                             3      6      2                             %
%                            o------o------o                              %           
%                           /             /                               %
%                        7 o             o 5                              %
%                         /             /                                 %
%                        o------o------o                                  %
%                        0      4      1                                  %
%                                                                         %
%-------------------------------------------------------------------------%
%                               Paraview HEXA 20                          %
%-------------------------------------------------------------------------%


%=========================================================================%
%                          INFOS MAILLAGE                                 %
%=========================================================================%

nb_node = length(nodes);                            % | Nb nodes
[nb_elem,taille_elem] = size(elems_para);           % | Nb elem and Nb node per elems
vector_taille_elem = taille_elem*ones(nb_elem,1);   % | Vecteur Nb node per elems
vtk_elem = [vector_taille_elem elems_para];         % | For CELLS in VTK
cell_types = 23*ones(nb_elem,1);                    % | For CELL_TYPES in VTK

%=========================================================================%
%                        FILE GENERATION .VTK                             %
%=========================================================================%

vtkfile = fopen([filename,'.vtk'],'w');
fprintf(vtkfile,'%s\n','# vtk DataFile Version 3.0');
fprintf(vtkfile,'%s\n','Unstructured Grid Example');
fprintf(vtkfile,'%s\n','ASCII');
fprintf(vtkfile,'%s\n','DATASET UNSTRUCTURED_GRID');
fprintf(vtkfile,'%s\n',['POINTS ',num2str(nb_node),' float']);
fprintf(vtkfile,'%i\t %i\t %i\n',nodes');
fprintf(vtkfile,'%s\n',['CELLS ', num2str(nb_elem),' ',num2str(nb_elem*9)]);
fprintf(vtkfile,'%i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\n',vtk_elem');
fprintf(vtkfile,'%s\n',['CELL_TYPES ', num2str(nb_elem)]);
fprintf(vtkfile,'%i\n',cell_types');
fprintf(vtkfile,'%s\n',' ');
fprintf(vtkfile,'%s\n',['POINT_DATA ',num2str(nb_node)]);
fprintf(vtkfile,'%s\n','SCALARS pres_adim float 1');
fprintf(vtkfile,'%s\n','LOOKUP_TABLE default');
fprintf(vtkfile,'%i\n',P');
fprintf(vtkfile,'%s\n',' ');
fprintf(vtkfile,'%s\n','VECTORS displacement float');
fprintf(vtkfile,'%i\t %i\t %i\n',U');
fclose(vtkfile);
disp(['File',filename,'.VTK OK'])

end

