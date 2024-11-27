function B = strainDispMatrix3d(nn,dRdx)
% Compute the 3D strain-displacement matrix
% using the order [epsilon_xx epsilon_yy epsilon_zz epsilon_xy epsilon_yz
% epsilon_zx] with U =[u_x^1
%                      u_x^n 
%                      u_y^1 
%                      u_y^n 
%                      u_z^1 u_z^n]
% WROTE BY Vinh Phu Nguyen
% Delft University of Technology

B(1,1:nn)         = dRdx(:,1)';
B(2,nn+1:2*nn)    = dRdx(:,2)';
B(3,2*nn+1:3*nn)  = dRdx(:,3)';

B(4,1:nn)         = dRdx(:,2)';
B(4,nn+1:2*nn)    = dRdx(:,1)';

B(5,2*nn+1:3*nn)  = dRdx(:,2)';
B(5,nn+1:2*nn)    = dRdx(:,3)';

B(6,1:nn)         = dRdx(:,3)';
B(6,2*nn+1:3*nn)  = dRdx(:,1)';
