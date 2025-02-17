function [IGA] = Car_Cavity()
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

p=2;
t = 0.035;
q = 2;
r = 2;
wknot = [0 0 0 1 1 1];
z0    = 0;
zL    = 1.2; 
zL2   = zL/2;

uKnot_1 = [0 0 0 1 1 1];
B_1= [0.1   -0.15;
      0.1   -0.075;
      0.1   -0.0];
weights_1 = ones(length(B_1),1);
figure
% plotNurbs(uKnot_1,p,weights_1,300,B_1)

uKnot_2 = [0 0 0 1 1 1];
B_2     = [0.1   -0.0;
           0.1   0.075;
           0.1   0.15
           ];
weights_2 = ones(length(B_2),1);
% plotNurbs(uKnot_2,p,weights_2,300,B_2)

uKnot_3 = [0 0 0 1 1 1];
B_3 = [0.1   0.15;
       0.1   0.225;
       0.1   0.3
      ];
weights_3 = ones(length(B_3),1);
% plotNurbs(uKnot_3,p,weights_3,300,B_3)

uKnot_4 = [0 0 0 1 1 1];
B_4 = [0.1   0.3;
       0.25  0.3; 
       0.4   0.3;
      ];
weights_4 = ones(length(B_4),1);
% plotNurbs(uKnot_4,p,weights_4,300,B_4)

uKnot_5 = [0 0 0 1 1 1];
B_5 = [0.4   0.3;
       0.85  0.3; 
       1   0.55;
       ];
weights_5 = ones(length(B_5),1);
% plotNurbs(uKnot_5,p,weights_5,300,B_5)

uKnot_6 = [0 0 0 1 1 1];
B_6 = [1   0.55;
       1.2  1;
       2.8  1
       ];
weights_6 = ones(length(B_6),1);
% plotNurbs(uKnot_6,p,weights_6,300,B_6)

uKnot_7 = [0 0 0 1 1 1];
B_7 = [2.8  1;
       3.5    1;
       3.5    0.5;
       ];
weights_7 = ones(length(B_7),1);
% plotNurbs(uKnot_7,p,weights_7,300,B_7)

uKnot_8 = [0 0 0 1 1 1];
B_8 = [3.5    0.5;
       3.5   -0.15;    
       2.4  -0.150;
       ];
weights_8 = ones(length(B_8),1);
% plotNurbs(uKnot_8,p,weights_8,300,B_8)

uKnot_9 = [0 0 0 1 1 1];
B_9 = [2.4  -0.150;
       1.4   -0.150; 
       0.4   -0.150;   
      ];
weights_9 = ones(length(B_9),1);
% plotNurbs(uKnot_9,p,weights_9,300,B_9)

uKnot_10 = [0 0 0 1 1 1];
B_10 = [0.1   -0.150;    
        0.25 -0.150; 
        0.4   -0.150;
      ];
weights_10 = ones(length(B_10),1);
% figure
% plotNurbs(uKnot_10,p,weights_10,300,B_10)

uKnot_11 = [0 0 0 1 1 1];
B_11 = [0.4   0.3;
        0.4   0.225; 
        0.4   0.15       
      ];
weights_11 = ones(length(B_11),1);
% plotNurbs(uKnot_11,p,weights_11,300,B_11)

uKnot_12 = [0 0 0 1 1 1];
B_12 = [0.4   0.;
        0.4   0.075; 
        0.4   0.15       
      ];
weights_12 = ones(length(B_12),1);
% plotNurbs(uKnot_12,p,weights_12,300,B_12)

uKnot_13 = [0 0 0 1 1 1];
B_13 = [0.4   -0.15;
        0.4   -0.075; 
        0.4   -0.0       
      ];
weights_13 = ones(length(B_13),1);
% plotNurbs(uKnot_13,p,weights_13,300,B_13)

uKnot_14 = [0 0 0 1 1 1];
B_14 = [0.1   0.15;
        0.25   0.15; 
        0.4   0.15       
      ];
weights_14 = ones(length(B_14),1);
% plotNurbs(uKnot_14,p,weights_14,300,B_14)

uKnot_15 = [0 0 0 1 1 1];
B_15 = [0.1   0.0;
        0.25   0.0; 
        0.4   0.0       
      ];
weights_15 = ones(length(B_15),1);
% plotNurbs(uKnot_15,p,weights_15,300,B_15)

uKnot_16 = [0 0 0 1 1 1];
B_16 = [1.6  0.4;
        1.3  0.475; 
        1    0.55       
      ];
weights_16 = ones(length(B_16),1);
% plotNurbs(uKnot_16,p,weights_16,300,B_16)

uKnot_17 = [0 0 0 1 1 1];
B_17 = [2.9    0.7;
        2.85   0.85; 
        2.8    1    
      ];
weights_17 = ones(length(B_17),1);
% plotNurbs(uKnot_17,p,weights_17,300,B_17)

uKnot_18 = [0 0 0 1 1 1];
B_18 = [3.5    0.5;
        3.3    0.5; 
        3.1    0.5       
      ];
weights_18 = ones(length(B_18),1);
% plotNurbs(uKnot_18,p,weights_18,300,B_18)

uKnot_19 = [0 0 0 1 1 1];
B_19 = [1.6    0.4;
        2.25  0.55; 
        2.9    0.7;           
      ];
weights_19 = ones(length(B_19),1);
% plotNurbs(uKnot_19,p,weights_19,300,B_19)

uKnot_20 = [0 0 0 1 1 1];
B_20 = [2.9    0.7;
        3      0.6; 
        3.1    0.5       
      ];
weights_20 = ones(length(B_20),1);
% figure
% plotNurbs(uKnot_20,p,weights_20,300,B_20)

uKnot_21 = [0 0 0 1 1 1];
B_21 = [0.4   0.15;
        1.3   0.15; 
        1.6   0.4       
      ];
weights_21 = ones(length(B_21),1);
% figure
% plotNurbs(uKnot_21,p,weights_21,300,B_21)

uKnot_22 = [0 0 0 1 1 1];
B_22 = [2    0.2;
        1.8   0.3; 
        1.6   0.4       
      ];
weights_22 = ones(length(B_22),1);
% plotNurbs(uKnot_22,p,weights_22,300,B_22)

uKnot_23 = [0 0 0 1 1 1];
B_23 = [2    0.2;
        2.55   0.35; 
        3.1   0.5       
      ];
weights_23 = ones(length(B_23),1);
% figure
% plotNurbs(uKnot_23,p,weights_23,300,B_23)

uKnot_24 = [0 0 0 1 1 1];
B_24 = [2.4    -0.15;
        2.2  0.025; 
        2    0.2       
      ];
weights_24 = ones(length(B_24),1);
% plotNurbs(uKnot_24,p,weights_24,300,B_24)

uKnot_25 = [0 0 0 1 1 1];
B_25 = [0.4  0.0;
        1.4    0.0; 
        2    0.2       
        ];
weights_25 = ones(length(B_25),1);
% plotNurbs(uKnot_25,p,weights_25,300,B_25)

% figure
% for i=1:19
%     uKnot_var_name = sprintf('uKnot_%d', i);
%     weights_var_name = sprintf('weights_%d', i);
%     B_var_name = sprintf('B_%d', i);
%     current_uKnot = eval(uKnot_var_name);
%     current_weights = eval(weights_var_name);
%     current_B = eval(B_var_name);
%     plotNurbs(current_uKnot,p,current_weights,300,current_B)
%     axis equal
%     hold on
%     % keyboard
% end
% hold on
% for i=20:25
%     uKnot_var_name = sprintf('uKnot_%d', i);
%     weights_var_name = sprintf('weights_%d', i);
%     B_var_name = sprintf('B_%d', i);
%     current_uKnot = eval(uKnot_var_name);
%     current_weights = eval(weights_var_name);
%     current_B = eval(B_var_name);
%     plotNurbs(current_uKnot,p,current_weights,300,current_B)
%     axis equal
%     hold on
%     % keyboard
% end

% % GENERATE SURFACES

%--FIRST SURFACE

uknot_1_s = [0 0 0 1 1 1]; 
vknot_1_s = [0 0 0 1 1 1];

B_1_s   = [B_3(1,:)  ;  B_14(2,:)                                      ; B_14(3,:);
           B_3(2,:)  ; [(B_3(2,1)+B_11(2,1))/2,(B_3(2,2)+B_11(2,2))/2]  ; B_11(2,:); 
           B_3(3,:)  ;  B_4(2,:)                                      ; B_4(3,:)];

weights_1_s = ones(size(B_1_s,1),1);
% figure
% plot(B_1_s(:,1),B_1_s(:,2),'ob','markersize',3,'markerface','b')
% for i=1:size(B_1_s,1)
%     text(B_1_s(i,1),B_1_s(i,2),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Color','b')
% end
%--1ST VOLUME

uknot_1_v = uknot_1_s; 
vknot_1_v = vknot_1_s;
wknot_1_v = wknot; 
noPtsX_1_v = length(uknot_1_v)-p-1;
noPtsY_1_v = length(vknot_1_v)-q-1;
noPtsZ_1_v = length(wknot_1_v)-r-1;
B_1_v = [B_1_s,z0*ones(size(B_1_s,1),1);
         B_1_s,zL2*ones(size(B_1_s,1),1);
         B_1_s,zL*ones(size(B_1_s,1),1);];

weights_1_v = ones(size(B_1_v,1),1);

% ElemMesh_vtk(B_1_v, weights_1_v,uknot_1_v,vknot_1_v,wknot_1_v,noPtsX_1_v,noPtsY_1_v,noPtsZ_1_v,p,q,r,'CAR_CAV/ELEMENT1')
% NurbsSurface_vtk(noPtsX_1_v,noPtsY_1_v,noPtsZ_1_v,uknot_1_v,vknot_1_v,wknot_1_v,B_1_v,weights_1_v,p,q,r,'CAR_CAV/NURBS1')

%--SECOND SURFACE

uknot_2_s = [0 0 0 1 1 1]; 
vknot_2_s = [0 0 0 1 1 1];

B_2_s   = [ B_2(1,:) ;  B_15(2,:)                                     ; B_15(3,:);
            B_2(2,:) ; [(B_2(2,1)+B_12(2,1))/2,(B_2(2,2)+B_12(2,2))/2] ; B_12(2,:); 
            B_2(3,:) ;  B_14(2,:)                                   ; B_12(3,:)];

weights_2_s = ones(size(B_2_s,1),1);
% figure
% plot(B_2_s(:,1),B_2_s(:,2),'ob','markersize',3,'markerface','b')
% for i=1:size(B_2_s,1)
%     text(B_2_s(i,1),B_2_s(i,2),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Color','b')
% end
%--2ND VOLUME

uknot_2_v = uknot_2_s; 
vknot_2_v = vknot_2_s;
wknot_2_v = wknot; 
noPtsX_2_v = length(uknot_2_v)-p-1;
noPtsY_2_v = length(vknot_2_v)-q-1;
noPtsZ_2_v = length(wknot_2_v)-r-1;
% B_2_v = [B_2_s,z0*ones(size(B_2_s,1),1);
%          B_2_s,zL2*ones(size(B_2_s,1),1);
%          B_2_s,zL*ones(size(B_2_s,1),1);];
B_2_v = [B_2_s,z0*ones(size(B_2_s,1),1);
         B_2_s,zL2*ones(size(B_2_s,1),1);
         B_2_s,zL*ones(size(B_2_s,1),1);];

% keyboard
weights_2_v = ones(size(B_2_v,1),1);
% 
% ElemMesh_vtk(B_2_v, weights_2_v,uknot_2_v,vknot_2_v,wknot_2_v,noPtsX_2_v,noPtsY_2_v,noPtsZ_2_v,p,q,r,'CAR_CAV/ELEMENT2')
% NurbsSurface_vtk(noPtsX_2_v,noPtsY_2_v,noPtsZ_2_v,uknot_2_v,vknot_2_v,wknot_2_v,B_2_v,weights_2_v,p,q,r,'CAR_CAV/NURBS2')

%--THIRD SURFACE

uknot_3_s = [0 0 0 1 1 1]; 
vknot_3_s = [0 0 0 1 1 1];

B_3_s = [ B_1(1,:) ; B_10(2,:)                                     ; B_10(3,:);
          B_1(2,:)  ; [(B_1(2,1)+B_13(2,1))/2,(B_1(2,2)+B_13(2,2))/2] ; B_13(2,:); 
          B_1(3,:)  ;  B_15(2,:)                                     ; B_15(3,:)];

weights_3_s = ones(size(B_3_s,1),1);
% figure
% plot(B_3_s(:,1),B_3_s(:,2),'ob','markersize',3,'markerface','b')
% for i=1:size(B_3_s,1)
%     text(B_3_s(i,1),B_3_s(i,2),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Color','b')
% end
%--3RD VOLUME

uknot_3_v = uknot_3_s; 
vknot_3_v = vknot_3_s;
wknot_3_v = wknot; 
noPtsX_3_v = length(uknot_3_v)-p-1;
noPtsY_3_v = length(vknot_3_v)-q-1;
noPtsZ_3_v = length(wknot_3_v)-r-1;
B_3_v = [B_3_s,z0*ones(size(B_3_s,1),1);
         B_3_s,zL2*ones(size(B_3_s,1),1);
         B_3_s,zL*ones(size(B_3_s,1),1);];

weights_3_v = ones(size(B_3_v,1),1);

% ElemMesh_vtk(B_3_v, weights_3_v,uknot_3_v,vknot_3_v,wknot_3_v,noPtsX_3_v,noPtsY_3_v,noPtsZ_3_v,p,q,r,'CAR_CAV/ELEMENT3')
% NurbsSurface_vtk(noPtsX_3_v,noPtsY_3_v,noPtsZ_3_v,uknot_3_v,vknot_3_v,wknot_3_v,B_3_v,weights_3_v,p,q,r,'CAR_CAV/NURBS3')

%--4TH SURFACE

uknot_4_s = [0 0 0 1 1 1]; 
vknot_4_s = [0 0 0 1 1 1];

B_4_s = [ B_5(3,:) ;    B_5(2,:)                                         ; B_5(1,:);
          B_16(2,:)  ; [(B_16(2,1)+B_11(2,1))/2+0.1,(B_16(2,2)+B_11(2,2))/2-0.1] ; B_11(2,:); 
          B_21(3,:)  ;   B_21(2,:)                                       ; B_21(1,:)];

weights_4_s = ones(size(B_4_s,1),1);
% figure
% plot(B_4_s(:,1),B_4_s(:,2),'ob','markersize',3,'markerface','b')
% for i=1:size(B_4_s,1)
%     text(B_4_s(i,1),B_4_s(i,2),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Color','b')
% end
%--4TH VOLUME

uknot_4_v = uknot_4_s; 
vknot_4_v = vknot_4_s;
wknot_4_v = wknot; 
noPtsX_4_v = length(uknot_4_v)-p-1;
noPtsY_4_v = length(vknot_4_v)-q-1;
noPtsZ_4_v = length(wknot_4_v)-r-1;
B_4_v = [B_4_s,z0*ones(size(B_4_s,1),1);
         B_4_s,zL2*ones(size(B_4_s,1),1);
         B_4_s,zL*ones(size(B_4_s,1),1);];

weights_4_v = ones(size(B_4_v,1),1);

% ElemMesh_vtk(B_4_v, weights_4_v,uknot_4_v,vknot_4_v,wknot_4_v,noPtsX_4_v,noPtsY_4_v,noPtsZ_4_v,p,q,r,'CAR_CAV/ELEMENT4')
% NurbsSurface_vtk(noPtsX_4_v,noPtsY_4_v,noPtsZ_4_v,uknot_4_v,vknot_4_v,wknot_4_v,B_4_v,weights_4_v,p,q,r,'CAR_CAV/NURBS4')

%--5TH SURFACE

uknot_5_s = [0 0 0 1 1 1]; 
vknot_5_s = [0 0 0 1 1 1];

B_5_s = [ B_21(3,:) ;    B_21(2,:)                                         ; B_21(1,:);
          B_22(2,:)  ; [(B_12(2,1)+B_22(2,1))/2+0.25,(B_12(2,2)+B_22(2,2))/2-0.1] ; B_12(2,:); 
          B_25(3,:)  ;   B_25(2,:)                                       ; B_25(1,:)];

weights_5_s = ones(size(B_5_s,1),1);
% figure
% plot(B_5_s(:,1),B_5_s(:,2),'ob','markersize',3,'markerface','b')
% for i=1:size(B_5_s,1)
%     text(B_5_s(i,1),B_5_s(i,2),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Color','b')
% end
%--5TH VOLUME

uknot_5_v = uknot_5_s; 
vknot_5_v = vknot_5_s;
wknot_5_v = wknot; 
noPtsX_5_v = length(uknot_5_v)-p-1;
noPtsY_5_v = length(vknot_5_v)-q-1;
noPtsZ_5_v = length(wknot_5_v)-r-1;
B_5_v = [B_5_s,z0*ones(size(B_5_s,1),1);
         B_5_s,zL2*ones(size(B_5_s,1),1);
         B_5_s,zL*ones(size(B_5_s,1),1);];

weights_5_v = ones(size(B_5_v,1),1);

% ElemMesh_vtk(B_5_v, weights_5_v,uknot_5_v,vknot_5_v,wknot_5_v,noPtsX_5_v,noPtsY_5_v,noPtsZ_5_v,p,q,r,'CAR_CAV/ELEMENT5')
% NurbsSurface_vtk(noPtsX_5_v,noPtsY_5_v,noPtsZ_5_v,uknot_5_v,vknot_5_v,wknot_5_v,B_5_v,weights_5_v,p,q,r,'CAR_CAV/NURBS5')

%--6TH SURFACE

uknot_6_s = [0 0 0 1 1 1]; 
vknot_6_s = [0 0 0 1 1 1];

B_6_s = [ B_25(3,:) ;    B_25(2,:)                                         ; B_25(1,:);
          B_24(2,:)  ; [(B_13(2,1)+B_24(2,1))/2+0.25,(B_13(2,2)+B_24(2,2))/2-0.01] ; B_13(2,:); 
          B_9(1,:)  ;   B_9(2,:)                                       ; B_9(3,:)];

weights_6_s = ones(size(B_6_s,1),1);
% figure
% plot(B_6_s(:,1),B_6_s(:,2),'ob','markersize',3,'markerface','b')
% for i=1:size(B_6_s,1)
%     text(B_6_s(i,1),B_6_s(i,2),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Color','b')
% end
%--6TH VOLUME

uknot_6_v = uknot_6_s; 
vknot_6_v = vknot_6_s;
wknot_6_v = wknot; 
noPtsX_6_v = length(uknot_6_v)-p-1;
noPtsY_6_v = length(vknot_6_v)-q-1;
noPtsZ_6_v = length(wknot_6_v)-r-1;
B_6_v = [B_6_s,z0*ones(size(B_6_s,1),1);
         B_6_s,zL2*ones(size(B_6_s,1),1);
         B_6_s,zL*ones(size(B_6_s,1),1);];

weights_6_v = ones(size(B_6_v,1),1);

% ElemMesh_vtk(B_6_v, weights_6_v,uknot_6_v,vknot_6_v,wknot_6_v,noPtsX_6_v,noPtsY_6_v,noPtsZ_6_v,p,q,r,'CAR_CAV/ELEMENT6')
% NurbsSurface_vtk(noPtsX_6_v,noPtsY_6_v,noPtsZ_6_v,uknot_6_v,vknot_6_v,wknot_6_v,B_6_v,weights_6_v,p,q,r,'CAR_CAV/NURBS6')

%--7TH SURFACE

uknot_7_s = [0 0 0 1 1 1]; 
vknot_7_s = [0 0 0 1 1 1];

B_7_s = [ B_19(1,:) ;    B_19(2,:)                                         ; B_19(3,:);
          B_16(2,:)  ; [(B_16(2,1)+B_17(2,1))/2,(B_16(2,2)+B_17(2,2))/2] ; B_17(2,:); 
          B_6(1,:)  ;   B_6(2,:)                                       ; B_6(3,:)];

weights_7_s = ones(size(B_7_s,1),1);
% figure
% plot(B_7_s(:,1),B_7_s(:,2),'ob','markersize',3,'markerface','b')
% for i=1:size(B_7_s,1)
%     text(B_7_s(i,1),B_7_s(i,2),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Color','b')
% end
%--7TH VOLUME

uknot_7_v = uknot_7_s; 
vknot_7_v = vknot_7_s;
wknot_7_v = wknot; 
noPtsX_7_v = length(uknot_7_v)-p-1;
noPtsY_7_v = length(vknot_7_v)-q-1;
noPtsZ_7_v = length(wknot_7_v)-r-1;
B_7_v = [B_7_s,z0*ones(size(B_7_s,1),1);
         B_7_s,zL2*ones(size(B_7_s,1),1);
         B_7_s,zL*ones(size(B_7_s,1),1);];

weights_7_v = ones(size(B_7_v,1),1);

% ElemMesh_vtk(B_7_v, weights_7_v,uknot_7_v,vknot_7_v,wknot_7_v,noPtsX_7_v,noPtsY_7_v,noPtsZ_7_v,p,q,r,'CAR_CAV/ELEMENT7')
% NurbsSurface_vtk(noPtsX_7_v,noPtsY_7_v,noPtsZ_7_v,uknot_7_v,vknot_7_v,wknot_7_v,B_7_v,weights_7_v,p,q,r,'CAR_CAV/NURBS7')

%--8TH SURFACE

uknot_8_s = [0 0 0 1 1 1]; 
vknot_8_s = [0 0 0 1 1 1];

B_8_s = [ B_23(1,:) ;    B_23(2,:)                                         ; B_23(3,:);
          B_22(2,:)  ; [(B_22(2,1)+B_20(2,1))/2,(B_22(2,2)+B_20(2,2))/2] ; B_20(2,:); 
          B_19(1,:)  ;   B_19(2,:)                                       ; B_19(3,:)];

weights_8_s = ones(size(B_8_s,1),1);
% figure
% plot(B_8_s(:,1),B_8_s(:,2),'ob','markersize',3,'markerface','b')
% for i=1:size(B_8_s,1)
%     text(B_8_s(i,1),B_8_s(i,2),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Color','b')
% end
%--8TH VOLUME

uknot_8_v = uknot_8_s; 
vknot_8_v = vknot_8_s;
wknot_8_v = wknot; 
noPtsX_8_v = length(uknot_8_v)-p-1;
noPtsY_8_v = length(vknot_8_v)-q-1;
noPtsZ_8_v = length(wknot_8_v)-r-1;
B_8_v = [B_8_s,z0*ones(size(B_8_s,1),1);
         B_8_s,zL2*ones(size(B_8_s,1),1);
         B_8_s,zL*ones(size(B_8_s,1),1);];

weights_8_v = ones(size(B_8_v,1),1);

% ElemMesh_vtk(B_8_v, weights_8_v,uknot_8_v,vknot_8_v,wknot_8_v,noPtsX_8_v,noPtsY_8_v,noPtsZ_8_v,p,q,r,'CAR_CAV/ELEMENT8')
% NurbsSurface_vtk(noPtsX_8_v,noPtsY_8_v,noPtsZ_8_v,uknot_8_v,vknot_8_v,wknot_8_v,B_8_v,weights_8_v,p,q,r,'CAR_CAV/NURBS8')

%--9TH SURFACE

uknot_9_s = [0 0 0 1 1 1]; 
vknot_9_s = [0 0 0 1 1 1];

B_9_s = [ B_20(1,:) ;    B_20(2,:)                                         ; B_20(3,:);
          B_17(2,:)  ; [(B_18(2,1)+B_17(2,1))/2,(B_17(2,2)+B_18(2,2))/2] ; B_18(2,:); 
          B_7(1,:)  ;   B_7(2,:)                                       ; B_7(3,:)];

weights_9_s = ones(size(B_9_s,1),1);
% figure
% plot(B_9_s(:,1),B_9_s(:,2),'ob','markersize',3,'markerface','b')
% for i=1:size(B_9_s,1)
%     text(B_9_s(i,1),B_9_s(i,2),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Color','b')
% end
%--9TH VOLUME

uknot_9_v = uknot_9_s; 
vknot_9_v = vknot_9_s;
wknot_9_v = wknot; 
noPtsX_9_v = length(uknot_9_v)-p-1;
noPtsY_9_v = length(vknot_9_v)-q-1;
noPtsZ_9_v = length(wknot_9_v)-r-1;
B_9_v = [B_9_s,z0*ones(size(B_9_s,1),1);
         B_9_s,zL2*ones(size(B_9_s,1),1);
         B_9_s,zL*ones(size(B_9_s,1),1);];

weights_9_v = ones(size(B_9_v,1),1);

% ElemMesh_vtk(B_9_v, weights_9_v,uknot_9_v,vknot_9_v,wknot_9_v,noPtsX_9_v,noPtsY_9_v,noPtsZ_9_v,p,q,r,'CAR_CAV/ELEMENT9')
% NurbsSurface_vtk(noPtsX_9_v,noPtsY_9_v,noPtsZ_9_v,uknot_9_v,vknot_9_v,wknot_9_v,B_9_v,weights_9_v,p,q,r,'CAR_CAV/NURBS9')

%--10TH SURFACE

uknot_10_s = [0 0 0 1 1 1]; 
vknot_10_s = [0 0 0 1 1 1];

B_10_s = [B_23(3,:) ;    B_23(2,:)                                      ; B_23(1,:);
          B_18(2,:)  ; [(B_18(2,1)+B_24(2,1))/2,(B_24(2,2)+B_18(2,2))/2] ; B_24(2,:); 
          B_8(1,:)  ;   B_8(2,:)                                         ; B_8(3,:)];

weights_10_s = ones(size(B_10_s,1),1);
% figure
% plot(B_10_s(:,1),B_10_s(:,2),'ob','markersize',3,'markerface','b')
% for i=1:size(B_10_s,1)
%     text(B_10_s(i,1),B_10_s(i,2),num2str(i),'FontSize',10,'FontWeight','bold','HorizontalAlignment','left','Color','b')
% end

%--10TH VOLUME

uknot_10_v = uknot_10_s; 
vknot_10_v = vknot_10_s;
wknot_10_v = wknot; 
noPtsX_10_v = length(uknot_10_v)-p-1;
noPtsY_10_v = length(vknot_10_v)-q-1;
noPtsZ_10_v = length(wknot_10_v)-r-1;
B_10_v = [B_10_s,z0*ones(size(B_10_s,1),1);
         B_10_s,zL2*ones(size(B_10_s,1),1);
         B_10_s,zL*ones(size(B_10_s,1),1);];

weights_10_v = ones(size(B_10_v,1),1);

% ElemMesh_vtk(B_10_v, weights_10_v,uknot_10_v,vknot_10_v,wknot_10_v,noPtsX_10_v,noPtsY_10_v,noPtsZ_10_v,p,q,r,'CAR_CAV/ELEMENT10')
% NurbsSurface_vtk(noPtsX_10_v,noPtsY_10_v,noPtsZ_10_v,uknot_10_v,vknot_10_v,wknot_10_v,B_10_v,weights_10_v,p,q,r,'CAR_CAV/NURBS10')

IGA = [];
for i = 1:10
    IGA(i).controlPts       = eval(['B_' num2str(i) '_v']);
    IGA(i).weights = eval(['weights_' num2str(i) '_v']);
    IGA(i).uknot       = eval(['uknot_' num2str(i) '_v']);
    IGA(i).vknot       = eval(['vknot_' num2str(i) '_v']);
    IGA(i).wknot       = eval(['wknot_' num2str(i) '_v']);
    IGA(i).p       = p;
    IGA(i).q       = q;
    IGA(i).r       = r;
    IGA(i).noPtsX   = eval(['noPtsX_' num2str(i) '_v']);
    IGA(i).noPtsY   = eval(['noPtsY_' num2str(i) '_v']);
    IGA(i).noPtsZ   = eval(['noPtsZ_' num2str(i) '_v']);
end
end

