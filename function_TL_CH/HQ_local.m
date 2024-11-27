function [Hl Ql]= HQ_local(X,Y,Z,c0,rhof)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% X: A vector containing the x-coordinates of the element's nodes.
% Y: A vector containing the y-coordinates of the element's nodes.
% Z: A vector containing the z-coordinates of the element's nodes.
% c0: A constant value used in the computation of the mass matrix.
% rhof: Density or another material property, used in the computation of the stiffness and mass matrices.

% OUTPUT:
% Hl: A sparse matrix representing the local stiffness matrix for the element, which is a 20×20 matrix.
% 
% Ql: A sparse matrix representing the local mass matrix for the element, also a 20×20 matrix.

Hl = sparse(zeros(20,20));
Ql = sparse(zeros(20,20));

[ir wp N dNeta dNnu dNte] = Point_intreg_3D(); 

for i=1:ir
    
    %----------------------------------------------------------------------
    % Jacobian of the deformation gradient
    %----------------------------------------------------------------------

    F(1,1)=dNeta(i,:)*X; F(1,2)=dNnu(i,:)*X; F(1,3)=dNte(i,:)*X;
    F(2,1)=dNeta(i,:)*Y; F(2,2)=dNnu(i,:)*Y; F(2,3)=dNte(i,:)*Y;
    F(3,1)=dNeta(i,:)*Z; F(3,2)=dNnu(i,:)*Z; F(3,3)=dNte(i,:)*Z;
    detF=det(F);

    if detF<10*eps
        disp('Jacobian determinant equal or less than zero!')
    end

    %----------------------------------------------------------------------
    % Derivatives of the shape functions
    %----------------------------------------------------------------------
    
    DD = [dNeta(i,:); dNnu(i,:); dNte(i,:)]';
    DN = DD/F;

    %----------------------------------------------------------------------
    % Discrete gradient operator and shape function matrix
    %----------------------------------------------------------------------
    
    He =[DN(1,1)  DN(2,1)  DN(3,1)  DN(4,1)  DN(5,1)  DN(6,1)  DN(7,1)  DN(8,1)  DN(9,1)  DN(10,1)  DN(11,1)  DN(12,1)  DN(13,1)  DN(14,1)  DN(15,1)  DN(16,1)  DN(17,1)  DN(18,1)  DN(19,1)  DN(20,1);
                DN(1,2)  DN(2,2)  DN(3,2)  DN(4,2)  DN(5,2)  DN(6,2)  DN(7,2)  DN(8,2)  DN(9,2)  DN(10,2)  DN(11,2)  DN(12,2)  DN(13,2)  DN(14,2)  DN(15,2)  DN(16,2)  DN(17,2)  DN(18,2)  DN(19,2)  DN(20,2);
                DN(1,3)  DN(2,3)  DN(3,3)  DN(4,3)  DN(5,3)  DN(6,3)  DN(7,3)  DN(8,3)  DN(9,3)  DN(10,3)  DN(11,3)  DN(12,3)  DN(13,3)  DN(14,3)  DN(15,3)  DN(16,3)  DN(17,3)  DN(18,3)  DN(19,3)  DN(20,3)];

    Qe =[N(i,1)  N(i,2)  N(i,3)  N(i,4)  N(i,5)  N(i,6)  N(i,7)  N(i,8)  N(i,9)  N(i,10)  N(i,11)  N(i,12)  N(i,13)  N(i,14)  N(i,15)  N(i,16)  N(i,17)  N(i,18)  N(i,19)  N(i,20) ];

    %----------------------------------------------------------------------
    % Stiffness and mass matrices
    %----------------------------------------------------------------------

    Hl= Hl + (1/(rhof)).*He'*He*detF*wp(i,1);
    Ql= Ql + (1/(rhof*c0*c0)).*Qe'*Qe*detF*wp(i,1);
    
end





