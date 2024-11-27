function [ output_args ] = matlab2VTK_net(nodes,elems,U,PID,filename)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT
% nodes: A matrix of node coordinates (N x 3), where each row corresponds to the (x, y, z) coordinates of a node.
% elems: A matrix defining the connectivity of the elements (M x N), where each row represents an element by listing the node indices that make it up.
% U: A vector or matrix containing the displacement data for the nodes.
% PID: A vector containing the property IDs associated with each element.
% filename: The name of the file (without extension) to which the VTK data will be written.
% OUTPUT:
% output_args: The function does not return any value explicitly. It primarily writes data to a .vtk file and displays a message indicating the successful creation of the file.

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

%=========================================================================%
%                          INFOS MAILLAGE                                 %
%=========================================================================%

nb_node = length(nodes);                            % | Nb nodes
[nb_elem,taille_elem] = size(elems_para);           % | Nb elem and Nb node per elems
vector_taille_elem = taille_elem*ones(nb_elem,1);   % | Vecteur Nb node per elems
vtk_elem = [vector_taille_elem elems_para];         % | For CELLS in VTK
cell_types = 3*ones(nb_elem,1);                    % | For CELL_TYPES in VTK

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
fprintf(vtkfile,'%s\n',['CELLS ', num2str(nb_elem),' ',num2str(nb_elem*3)]); %nnode per element +1
fprintf(vtkfile,'%i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\n',vtk_elem');
fprintf(vtkfile,'%s\n',' ');
fprintf(vtkfile,'%s\n',['CELL_TYPES ', num2str(nb_elem)]);
fprintf(vtkfile,'%i\n',cell_types');
fprintf(vtkfile,'%s\n',' ');
fprintf(vtkfile,'%s\n',['CELL_DATA ',num2str(nb_elem)]);
fprintf(vtkfile,'%s\n','SCALARS Pid float 1');
fprintf(vtkfile,'%s\n','LOOKUP_TABLE default');
fprintf(vtkfile,'%i\n',PID');
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

