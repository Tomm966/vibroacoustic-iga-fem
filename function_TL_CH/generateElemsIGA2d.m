function [nodesTot,elemsTot] = generateElemsIGA2d(controlPts_coupl, ...
                                             uKnot_coupl, ...
                                             vKnot_coupl, ...
                                             noPtsX_coupl, ...
                                             noPtsY_coupl, ...
                                             p,q,weights_coupl)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% controlPts_coupl: The control points for the NURBS surface.
% uKnot_coupl, vKnot_coupl: Knot vectors for the xi and eta directions.
% noPtsX_coupl, noPtsY_coupl: The number of control points in the xi and eta directions.
% p, q: The degrees of the NURBS in the respective directions.
% weights_coupl: Weights associated with the control points.
% OUTPUT:
% nodesTot: A matrix containing the physical node coordinates.
% elemsTot: A matrix detailing the connectivity of the elements.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Elements coupling matrix
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[element_coupl,elRangeU_coupl,...
 elConnU_coupl,elRangeV_coupl,...
 elConnV_coupl,noElems_coupl,index_coupl] = generateIGA2DMesh_CH(uKnot_coupl, ...
                                        vKnot_coupl, ...
                                        noPtsX_coupl,...
                                        noPtsY_coupl, ...
                                        p,q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


elemsTot = [];
nodesTot = [];
for e=1:noElems_coupl


   idu    = index_coupl(e,1);
   idv    = index_coupl(e,2);
   xiE    = elRangeU_coupl(idu,:); % [xi_i,xi_i+1]
   etaE   = elRangeV_coupl(idv,:); % [eta_j,eta_j+1]
    
   nnl = 10;
   nodes1  = [linspace(xiE(1),xiE(2),nnl)' , ones(nnl,1)*etaE(1)];
   nodes2  = [linspace(xiE(1),xiE(2),nnl)' , ones(nnl,1)*etaE(2)];
   nodes3  = [ones(nnl,1)*xiE(1) , linspace(etaE(1),etaE(2),nnl)'];
   nodes4  = [ones(nnl,1)*xiE(2) , linspace(etaE(1),etaE(2),nnl)'];
   nodes = [nodes1;nodes2;nodes3;nodes4];
   elems = zeros(4*(nnl-1),2); 
   for i = 1:4
        elemsi = [1:nnl-1 ; 2:nnl]';
        elems([1:nnl-1] + (i-1)*(nnl -1),:) = elemsi + (i-1)*nnl;
   end
   dim = 4;
   nodes_sur = zeros(size(nodes,1),dim);
   projcoord = nurb2proj(noPtsX_coupl*noPtsY_coupl,controlPts_coupl,weights_coupl);
   for j = 1:size(nodes,1)
      nodes_sur(j,:) = SurfacePoint(noPtsX_coupl-1,p, uKnot_coupl, ... 
                                    noPtsY_coupl-1,q, vKnot_coupl, ...
                                    projcoord,dim, ...
                                    nodes(j,1),nodes(j,2));
   end
   nodes_sur(:,[1,2,3]) = nodes_sur(:,[1,2,3])./nodes_sur(:,4); 
   elemsTot = [elemsTot ; elems + size(nodesTot,1)];
   nodesTot = [nodesTot ; nodes_sur];
end
end
