function [ output_args ] = matlab2VTK_Quad4(nodes,elems,U,P,filename)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% nodes: An N x 3 matrix containing the coordinates of each node in the mesh, where N is the number of nodes.
% elems: An M x 4 matrix defining the connectivity of each quadrilateral element, where each row specifies the indices of the four nodes that make up that element.
% U: A vector or matrix containing the displacement values for each node, reshaped to ensure it matches the format (N x 3).
% P: A vector containing scalar data (e.g., pressure) for each node, which can be visualized alongside the mesh.
% filename: The name of the output file (without extension) to which the VTK data will be written.
% OUTPUT:
% output_args: The function does not explicitly return any values; it primarily writes data to a .vtk file and confirms the successful creation of the file via a console message.

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
cell_types = 9*ones(nb_elem,1);                    % | For CELL_TYPES in VTK

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
fprintf(vtkfile,'%s\n',['CELLS ', num2str(nb_elem),' ',num2str(nb_elem*5)]);
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

