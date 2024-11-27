function [ output_args ] = matlab2VTK_points(nodes,elems,U,filename)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% nodes: A matrix of node coordinates, where each row corresponds to a (x, y, z) coordinate of a node (N x 3).
% elems: A matrix defining the connectivity of the elements, where each row specifies the nodes that make up an element. However, in this version of the function, elements are not used for point data.
% U: A vector or matrix containing the displacement data for the nodes (N x 3).
% filename: The name of the file (without extension) to which the VTK data will be written.
% OUTPUT:
% output_args: The function does not return any value explicitly; it primarily writes data to a .vtk file and confirms successful file creation via a message.

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
cell_types = 1*ones(nb_elem,1);                    % | For CELL_TYPES in VTK

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
fprintf(vtkfile,'%s\n',['CELLS ', num2str(nb_elem),' ',num2str(nb_elem*2)]); %nnode per element +1
fprintf(vtkfile,'%i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\t %i\n',vtk_elem');
fprintf(vtkfile,'%s\n',' ');
fprintf(vtkfile,'%s\n',['CELL_TYPES ', num2str(nb_elem)]);
fprintf(vtkfile,'%i\n',cell_types');
fprintf(vtkfile,'%s\n',' ');
% fprintf(vtkfile,'%s\n',['CELL_DATA ',num2str(nb_elem)]);
% fprintf(vtkfile,'%s\n','SCALARS Pid float 1');
% fprintf(vtkfile,'%s\n','LOOKUP_TABLE default');
% fprintf(vtkfile,'%i\n',PID');
% fprintf(vtkfile,'%s\n',' ');
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

