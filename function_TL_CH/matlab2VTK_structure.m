function [ output_args ] = matlab2VTK_structure(nodes,elems,U,filename)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% Inputs
% nodes: An n×3 matrix containing the x,y,z coordinates of the mesh nodes.
% elems: An m×N matrix where each row represents an element by listing the indices of its nodes.
% U: A vector or matrix representing the displacement of each node. It is reshaped into a n×3 matrix, where each row corresponds to the displacement in the x,y,z directions for each node.
% filename: The desired name of the output VTK file (without the .vtk extension).

% Outputs
% The function generates a .vtk file with the specified filename containing the mesh and displacement data.

nnode = length(nodes);
U = reshape(U,[3,nnode])';
%F = reshape(F,[3,nnode])';
%Fint = reshape(Fint,[3,nnode])';
%=========================================================================%
%                    VISUALIZATION WITH PARAVIEW                          %
%=========================================================================%

%-------------------------------------------------------------------------%
%         Matlab read data from 1 to N                                    %
%         Paraview read data from 0 to N-1                                %
%-------------------------------------------------------------------------%

elems_para = elems - ones(size(elems)); % Renumerotation paraview

%-------------------------------------------------------------------------%
%                               Nastran HEXA 20                           %
%-------------------------------------------------------------------------%
%                                                                         %
%                            7     18      6                              %
%                            o------o------o                              %
%                           /|            /|                              %
%                      19  o | 15     17 o |                              %
%                         /  o   16     /  o  14                          %
%                      4 o------o------o 5 |                              %
%                        |   |     10  |   |                              %
%                        | 3 o------o--|---o 2                            %           
%                     12 o  /          o  /                               %
%                        | o 11     13 | o  9                             %
%                        |/            |/                                 %
%                        o------o------o                                  %
%                        0      8      1                                  %
%                                                                         %
%-------------------------------------------------------------------------%
%                               Paraview HEXA 20                          %
%-------------------------------------------------------------------------%
%                                                                         %
%                            7     14      6                              %
%                            o------o------o                              %
%                           /|            /|                              %
%                      15  o | 19     13 o |                              %
%                         /  o   12     /  o  18                          %
%                      4 o------o------o 5 |                              %
%                        |   |     10  |   |                              %
%                        | 3 o------o--|---o 2                            %           
%                     16 o  /          o  /                               %
%                        | o 11     17 | o  9                             %
%                        |/            |/                                 %
%                        o------o------o                                  %
%                        0      8      1                                  %
%                                                                         %
%-------------------------------------------------------------------------%

stock_elem_1 = elems_para(:,(13:16));
stock_elem_2 = elems_para(:,(17:20));

% We have to change the numerotation for the visualization

elems_para(:,(13:16)) = stock_elem_2 ;
elems_para(:,(17:20)) = stock_elem_1 ;

%=========================================================================%
%                          INFOS MAILLAGE                                 %
%=========================================================================%

nb_node = length(nodes);                            % | Nb nodes
[nb_elem,taille_elem] = size(elems_para);           % | Nb elem and Nb node per elems
vector_taille_elem = taille_elem*ones(nb_elem,1);   % | Vecteur Nb node per elems
vtk_elem = [vector_taille_elem elems_para];         % | For CELLS in VTK
cell_types = 25*ones(nb_elem,1);                    % | For CELL_TYPES in VTK

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
fprintf(vtkfile,'%s\n',' ');
fprintf(vtkfile,'%s\n',['CELLS ', num2str(nb_elem),' ',num2str(nb_elem*21)]);
fprintf(vtkfile,'%i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\n',vtk_elem');
fprintf(vtkfile,'%s\n',' ');
fprintf(vtkfile,'%s\n',['CELL_TYPES ', num2str(nb_elem)]);
fprintf(vtkfile,'%i\n',cell_types');
fprintf(vtkfile,'%s\n',' ');
fprintf(vtkfile,'%s\n',['POINT_DATA ',num2str(nb_node)]);
% fprintf(vtkfile,'%s\n',['SCALARS Epsilon float ',num2str(6)]);
% fprintf(vtkfile,'%s\n','LOOKUP_TABLE default');
% fprintf(vtkfile,'%i\t %i\t %i\t %i\t %i\t %i\n',Eps(:,1:6)');
% fprintf(vtkfile,'%s\n',' ');
fprintf(vtkfile,'%s\n','VECTORS displacement float');
fprintf(vtkfile,'%i\t %i\t %i\n',U');

% fprintf(vtkfile,'%s\n',' ');
% fprintf(vtkfile,'%s\n','VECTORS Fext float');
% fprintf(vtkfile,'%i\t %i\t %i\n',F');
% fprintf(vtkfile,'%s\n',' ');
% fprintf(vtkfile,'%s\n','VECTORS Fint float');
% fprintf(vtkfile,'%i\t %i\t %i\n',Fint');
fclose(vtkfile);
disp(['File',filename,'.VTK OK'])

end

