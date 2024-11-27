function [Ce,Area]=C_local_courbe(X,Y,Z,rhof)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

%INPUT:
% X: Coordinate matrix of the element in the x-direction (vector).
% Y: Coordinate matrix of the element in the y-direction (vector).
% Z: Coordinate matrix of the element in the z-direction (vector).
% rhof: Density of the fluid (scalar).

% OUTPUT:
% Ce: Element coupling matrix (24x8 matrix).
% Area: Surface area of the element (scalar).


%=========================================================================%
%                     Gauss weight and Gauss point (cf. T.GMUR p 231)     %
%=========================================================================%

ir=3;
ngp=ir*ir;

if ir==1
    g1=0.0; w1=2.0;
    gp=[ g1 g1 ];  w=[ w1 w1 ];
elseif ir==2
    g1=0.577350269189626; w1=1;
    gp(:,1)=[-g1; g1;-g1; g1];  gp(:,2)=[-g1;-g1; g1; g1];
    w(:,1)=[ w1; w1; w1; w1];   w(:,2)=[ w1; w1; w1; w1];
elseif ir==3
    g1=0.774596669241483; g2=0.;
    w1=0.555555555555555; w2=0.888888888888888;
    gp(:,1)=[-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2)=[-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    w(:,1)=[ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2)=[ w1; w1; w1; w2; w2; w2; w1; w1; w1];
else
    disp('attention au nombre de points d''integration ');
    return
end
wp=w(:,1).*w(:,2);

xsi=gp(:,1);  eta=gp(:,2);  r2=ngp*2;

%=========================================================================%
%                             Shape functions                             %
%=========================================================================%

N(:,1)=(1-xsi).*(1-eta)/4;
N(:,2)=(1+xsi).*(1-eta)/4;
N(:,3)=(1+xsi).*(1+eta)/4;
N(:,4)=(1-xsi).*(1+eta)/4;
N(:,5)=(1-xsi.*xsi).*(1-eta)/2;
N(:,6)=(1+xsi).*(1-eta.*eta)/2;
N(:,7)=(1-xsi.*xsi).*(1+eta)/2;
N(:,8)=(1-xsi).*(1-eta.*eta)/2;

N(:,1)=N(:,1) - ((N(:,8)+N(:,5))/2);
N(:,2)=N(:,2) - ((N(:,5)+N(:,6))/2);
N(:,3)=N(:,3) - ((N(:,6)+N(:,7))/2);
N(:,4)=N(:,4) - ((N(:,7)+N(:,8))/2);

%=========================================================================%
%                       Shape functions /xsi dérivation                   %
%=========================================================================%

dNr(1:2:r2,1)=-(1-eta)/4;
dNr(1:2:r2,2)= (1-eta)/4;
dNr(1:2:r2,3)= (1+eta)/4;           
dNr(1:2:r2,4)=-(1+eta)/4;
dNr(1:2:r2,5)=-xsi.*(1-eta);
dNr(1:2:r2,6)=(1-eta.*eta)/2;
dNr(1:2:r2,7)=-xsi.*(1+eta);
dNr(1:2:r2,8)=-(1-eta.*eta)/2;

dNr(1:2:r2,1)= dNr(1:2:r2,1) - ((dNr(1:2:r2,8) + dNr(1:2:r2,5))/2) ;
dNr(1:2:r2,2)= dNr(1:2:r2,2) - ((dNr(1:2:r2,5) + dNr(1:2:r2,6))/2) ;
dNr(1:2:r2,3)= dNr(1:2:r2,3) - ((dNr(1:2:r2,6) + dNr(1:2:r2,7))/2) ;            
dNr(1:2:r2,4)= dNr(1:2:r2,4) - ((dNr(1:2:r2,7) + dNr(1:2:r2,8))/2) ;


%=========================================================================%
%                       Shape functions /eta dérivation                   %
%=========================================================================%


dNr(2:2:r2+1,1)=-(1-xsi)/4; 
dNr(2:2:r2+1,2)=-(1+xsi)/4;
dNr(2:2:r2+1,3)= (1+xsi)/4;
dNr(2:2:r2+1,4)= (1-xsi)/4;
dNr(2:2:r2+1,5)= -(1-xsi.*xsi)/2;
dNr(2:2:r2+1,6)= -(1+xsi).*eta;
dNr(2:2:r2+1,7)= (1-xsi.*xsi)/2;
dNr(2:2:r2+1,8)= (1-xsi).*eta;

dNr(2:2:r2+1,1)= dNr(2:2:r2+1,1) - ((dNr(2:2:r2+1,8) + dNr(2:2:r2+1,5))/2) ;
dNr(2:2:r2+1,2)= dNr(2:2:r2+1,2) - ((dNr(2:2:r2+1,5) + dNr(2:2:r2+1,6))/2) ;
dNr(2:2:r2+1,3)= dNr(2:2:r2+1,3) - ((dNr(2:2:r2+1,6) + dNr(2:2:r2+1,7))/2) ;            
dNr(2:2:r2+1,4)= dNr(2:2:r2+1,4) - ((dNr(2:2:r2+1,7) + dNr(2:2:r2+1,8))/2) ;


%=========================================================================%
%                       MATRIX WITH SHAPE FNUCTIONS  (W. Larbi p 175)     %
%=========================================================================%
%
%----------------- Displacement shape function matrix (3x24) -------------%
%                                                                         %
%           _                                _                            %
%     NU = |  N1 0  0  N2 0  0  ... N8 0  0   |                           %
%          |  0  N1 0  0  N2 0  ... 0  N8 0   |                           %
%          |_ 0  0  N1 0  0  N2 ... 0  0  N8 _|                           %                                      
%                                                                         %
%                                                                         %
%------------------- Pressure shape function matrix (1x8) ----------------%
%                                                                         %
%           _                  _                                          %
%     NP = |_  N1 N2  ...  N8  _|                                         %
%                                                                         %
%                                                                         %
%------------------- Elementary coupling matrix (24x8) -------------------%
%                                                                         %
%                                                                         %
%     Ce = NU'*n*NP*detF*wg                                               %
%           _          _     _    _     _                  _              %    
%        = |  N1 0  0   | * |  n1  | * |_  N1 N2  ...  N8  _|*detF*wg     %
%          |  0  N1 0   |   |  n2  |                                      %
%          |  0  0  N1  |   |_ n3 _|                                      %
%          |  .  .  .   |                                                 %
%          |  .  .  .   |                                                 %
%          |_ 0  0  N8 _|                                                 %
%                                                                         %
%                                                                         %
%=========================================================================%

NU=zeros(3,24);
NP=zeros(1,8) ;
n = zeros(3,1);

Ce=zeros(size(NU'*n*NP));
Area = 0;

%=========================================================================%
%                         JACOBIAN CALCULATION                            %
%=========================================================================%

JT=dNr*[X,Y,Z];

%=========================================================================%
%                             Ce calculation                              %
%=========================================================================%

    for i=1:ngp
        indx=[ 2*i-1; 2*i ];
        J = JT(indx,:)';
        detJxy = J(1,1)*J(2,2) - J(2,1)*J(1,2);
        detJyz = J(2,1)*J(3,2) - J(3,1)*J(2,2);
        detJzx = J(3,1)*J(1,2) - J(1,1)*J(3,2);    
        detJ = sqrt(detJxy^2 + detJyz^2 + detJzx^2);

        tr = J(:,1);
        ts = J(:,2);
        n = (1/norm(cross(tr,ts),2)).*cross(tr,ts);

        if detJ<0
            disp('probleme de jacobien !')
        end
        NP(1,1:8)     =N(i,:);
        NU(1,1:3:24-2)=N(i,:);
        NU(2,2:3:24-1)=N(i,:);
        NU(3,3:3:24)  =N(i,:);
        Area = Area + detJ*wp(i);
        Ce= Ce+ NU'*n*NP*detJ*wp(i);

    end

%=========================================================================%

end


