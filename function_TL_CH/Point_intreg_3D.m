function [ir wp N dNeta dNnu dNte] = Point_intreg_3D()
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% Input
% ir: Numero di punti di integrazione (intero; può essere 8 o 27).

% Output
% ir: Numero di punti di integrazione (intero).
% wp: Pesi di integrazione (matrice 27x1 per 27 punti, matrice 8x1 per 8 punti).
% N: Funzioni di forma (matrice 20xN, dove N è il numero di punti di integrazione).
% dNeta: Derivata delle funzioni di forma rispetto a ETA (matrice 20xN).
% dNnu: Derivata delle funzioni di forma rispetto a NI (matrice 20xN).
% dNte: Derivata delle funzioni di forma rispetto a XSI (matrice 20xN).

ir=27;

%--------------------------------------------------------------------------
% Integration points - CUB20
%--------------------------------------------------------------------------

if ir==8 % 2 points in each direction (order 3) %a = 1/sqrt(3);
    a = 5.77350269189626e-001;
    gp(:,1)=[-a; -a; -a; -a; a; a; a; a];  
    gp(:,2)=[-a; -a; a; a; -a; -a; a; a];  
    gp(:,3)=[-a; a; -a; a; -a; a; -a; a];  
    wp(:,1)=[1; 1; 1; 1; 1; 1; 1; 1];  
elseif ir==27 % 3 points in each direction (order 5)
    a = 7.74596669241483e-001; c1 = 5.55555555555556e-001; c2 = 8.88888888888889e-001; %a = sqrt(3/5); c1 = 5/9; c2 = 8/9; 
    gp(:,1)=[-a; -a; -a; -a; -a; -a; -a; -a; -a; 0; 0; 0; 0; 0; 0; 0; 0; 0; a; a; a; a; a; a; a; a; a];  
    gp(:,2)=[-a; -a; -a; 0; 0; 0; a; a; a; -a; -a; -a; 0; 0; 0; a; a; a; -a; -a; -a; 0; 0; 0; a; a; a]; 
    gp(:,3)=[-a; 0; a; -a; 0; a; -a; 0; a; -a; 0; a; -a; 0; a; -a; 0; a; -a; 0; a; -a; 0; a; -a; 0; a];
    wp(:,1)=[c1^3; c1^2*c2; c1^3; c1^2*c2; c1*c2^2; c1^2*c2;  c1^3; c1^2*c2; c1^3; c1^2*c2; c1*c2^2; c1^2*c2; c1*c2^2; c2^3; c1*c2^2; c1^2*c2; c1*c2^2; c2*c1^2; c1^3; c2*c1^2; c1^3; c1^2*c2; c1*c2^2; c1^2*c2; c1^3; c1^2*c2; c1^3];
else
    disp('Used number of integration points not implemented');
    return
end

eta=gp(:,1);  nu=gp(:,2); te = gp(:,3);

%--------------------------------------------------------------------------
% Fonctions de forme
%--------------------------------------------------------------------------

N(:,1)=(1/8)*(1-eta).*(1-nu).*(1-te).*(-2-eta-nu-te);   N(:,2)=(1/8)*(1+eta).*(1-nu).*(1-te).*(-2+eta-nu-te);
N(:,3)=(1/8)*(1+eta).*(1+nu).*(1-te).*(-2+eta+nu-te);   N(:,4)=(1/8)*(1-eta).*(1+nu).*(1-te).*(-2-eta+nu-te);
N(:,5)=(1/8)*(1-eta).*(1-nu).*(1+te).*(-2-eta-nu+te);   N(:,6)=(1/8)*(1+eta).*(1-nu).*(1+te).*(-2+eta-nu+te);
N(:,7)=(1/8)*(1+eta).*(1+nu).*(1+te).*(-2+eta+nu+te);   N(:,8)=(1/8)*(1-eta).*(1+nu).*(1+te).*(-2-eta+nu+te);
N(:,9) =0.25*(1-eta.^2).*(1-nu).*(1-te);                N(:,10)=0.25*(1+eta).*(1-nu.^2).*(1-te);
N(:,11)=0.25*(1-eta.^2).*(1+nu).*(1-te);                N(:,12)=0.25*(1-eta).*(1-nu.^2).*(1-te);
N(:,13)=0.25*(1-eta).*(1-nu).*(1-te.^2);                N(:,14)=0.25*(1+eta).*(1-nu).*(1-te.^2);
N(:,15)=0.25*(1+eta).*(1+nu).*(1-te.^2);                N(:,16)=0.25*(1-eta).*(1+nu).*(1-te.^2);
N(:,17)=0.25*(1-eta.^2).*(1-nu).*(1+te);                N(:,18)=0.25*(1+eta).*(1-nu.^2).*(1+te);
N(:,19)=0.25*(1-eta.^2).*(1+nu).*(1+te);                N(:,20)=0.25*(1-eta).*(1-nu.^2).*(1+te);

%--------------------------------------------------------------------------
% Dérivés de fonctions de forme
%--------------------------------------------------------------------------

dNeta(:,1)=(1/8)*(1-nu).*(1-te).*(1+2*eta+nu+te);       dNeta(:,2)=(1/8)*(1-nu).*(1-te).*(-1+2*eta-nu-te);
dNeta(:,3)=(1/8)*(1+nu).*(1-te).*(-1+2*eta+nu-te);      dNeta(:,4)=(1/8)*(1+nu).*(1-te).*(1+2*eta-nu+te);
dNeta(:,5)=(1/8)*(1-nu).*(1+te).*(1+2*eta+nu-te);       dNeta(:,6)=(1/8)*(1-nu).*(1+te).*(-1+2*eta-nu+te);
dNeta(:,7)=(1/8)*(1+nu).*(1+te).*(-1+2*eta+nu+te);      dNeta(:,8)=(1/8)*(1+nu).*(1+te).*(1+2*eta-nu-te);
dNeta(:,9)=0.25*(-2*eta).*(1-nu).*(1-te);               dNeta(:,10)=0.25*(1-nu.^2).*(1-te);
dNeta(:,11)=0.25*(-2*eta).*(1+nu).*(1-te);              dNeta(:,12)=-0.25*(1-nu.^2).*(1-te);
dNeta(:,13)=-0.25*(1-nu).*(1-te.^2);                    dNeta(:,14)=0.25*(1-nu).*(1-te.^2);
dNeta(:,15)=0.25*(1+nu).*(1-te.^2);                     dNeta(:,16)=-0.25*(1+nu).*(1-te.^2);
dNeta(:,17)=0.25*(-2*eta).*(1-nu).*(1+te);              dNeta(:,18)=0.25*(1-nu.^2).*(1+te);
dNeta(:,19)=0.25*(-2*eta).*(1+nu).*(1+te);              dNeta(:,20)=-0.25*(1-nu.^2).*(1+te);

dNnu(:,1)=(1/8)*(1-eta).*(1-te).*(1+eta+2*nu+te);       dNnu(:,2)=(1/8)*(1+eta).*(1-te).*(1-eta+2*nu+te);
dNnu(:,3)=(1/8)*(1+eta).*(1-te).*(-1+eta+2*nu-te);      dNnu(:,4)=(1/8)*(1-eta).*(1-te).*(-1-eta+2*nu-te);
dNnu(:,5)=(1/8)*(1-eta).*(1+te).*(1+eta+2*nu-te);       dNnu(:,6)=(1/8)*(1+eta).*(1+te).*(1-eta+2*nu-te);
dNnu(:,7)=(1/8)*(1+eta).*(1+te).*(-1+eta+2*nu+te);      dNnu(:,8)=(1/8)*(1-eta).*(1+te).*(-1-eta+2*nu+te);
dNnu(:,9)=-0.25*(1-eta.^2).*(1-te);                     dNnu(:,10)=0.25*(1+eta).*(-2*nu).*(1-te);
dNnu(:,11)=0.25*(1-eta.^2).*(1-te);                     dNnu(:,12)=0.25*(1-eta).*(-2*nu).*(1-te);
dNnu(:,13)=-0.25*(1-eta).*(1-te.^2);                    dNnu(:,14)=-0.25*(1+eta).*(1-te.^2);
dNnu(:,15)=0.25*(1+eta).*(1-te.^2);                     dNnu(:,16)=0.25*(1-eta).*(1-te.^2);
dNnu(:,17)=-0.25*(1-eta.^2).*(1+te);                    dNnu(:,18)=0.25*(1+eta).*(-2*nu).*(1+te);
dNnu(:,19)=0.25*(1-eta.^2).*(1+te);                     dNnu(:,20)=0.25*(1-eta).*(-2*nu).*(1+te);

dNte(:,1)=(1/8)*(1-eta).*(1-nu).*(1+eta+nu+2*te);       dNte(:,2)=(1/8)*(1+eta).*(1-nu).*(1-eta+nu+2*te);
dNte(:,3)=(1/8)*(1+eta).*(1+nu).*(1-eta-nu+2*te);       dNte(:,4)=(1/8)*(1-eta).*(1+nu).*(1+eta-nu+2*te);
dNte(:,5)=(1/8)*(1-eta).*(1-nu).*(-1-eta-nu+2*te);      dNte(:,6)=(1/8)*(1+eta).*(1-nu).*(-1+eta-nu+2*te);
dNte(:,7)=(1/8)*(1+eta).*(1+nu).*(-1+eta+nu+2*te);      dNte(:,8)=(1/8)*(1-eta).*(1+nu).*(-1-eta+nu+2*te);
dNte(:,9)=-0.25*(1-eta.^2).*(1-nu);                     dNte(:,10)=-0.25*(1+eta).*(1-nu.^2);
dNte(:,11)=-0.25*(1-eta.^2).*(1+nu);                    dNte(:,12)=-0.25*(1-eta).*(1-nu.^2);
dNte(:,13)=0.25*(1-eta).*(1-nu).*(-2*te);               dNte(:,14)=0.25*(1+eta).*(1-nu).*(-2*te);
dNte(:,15)=0.25*(1+eta).*(1+nu).*(-2*te);               dNte(:,16)=0.25*(1-eta).*(1+nu).*(-2*te);
dNte(:,17)=0.25*(1-eta.^2).*(1-nu);                     dNte(:,18)=0.25*(1+eta).*(1-nu.^2);
dNte(:,19)=0.25*(1-eta.^2).*(1+nu);                     dNte(:,20)=0.25*(1-eta).*(1-nu.^2);

