function [noDofs,Ka,Ma] = HQ_matrices_3D_parallel(NURBS,noGPs,c0)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% NURBS: A struct containing information about the NURBS (Non-Uniform Rational B-Splines) geometry, including knot vectors, control points, weights, and polynomial degrees.
% noGPs: The number of Gauss points used for numerical integration.
% c0: FLUID SPEED
% OUTPUT
% noDofs: The total number of degrees of freedom (DOFs).
% Ka: The assembled global stiffness matrix (sparse format).
% Ma: The assembled global mass matrix (sparse format)



uKnot = NURBS.uknot;
vKnot = NURBS.vknot;
wKnot = NURBS.wknot;
noPtsX = NURBS.noPtsX;
noPtsY = NURBS.noPtsY;
noPtsZ = NURBS.noPtsZ;
controlPts = NURBS.controlPts;
p = NURBS.p;
q = NURBS.q;
r = NURBS.r;
weights = NURBS.weights;

noCtrPts   = noPtsX * noPtsY* noPtsZ;
noDofs     = noCtrPts;
uniqueUKnots   = unique(uKnot);
uniqueVKnots   = unique(vKnot);
uniqueWKnots   = unique(wKnot);

noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.
noElemsW       = length(uniqueWKnots)-1; % # of elements zeta dir.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% chan =
%        1 2 3 4
%        5 6 7 8
%        9 10 11 12
%        13 14 15 16
% for a 4x2x2 control points


chan  = zeros(noPtsZ,noPtsY,noPtsX);

count = 1;

for i=1:noPtsZ
    for j=1:noPtsY
        for k=1:noPtsX
            chan(i,j,k) = count;
            count       = count + 1;
        end
    end
end


% determine our element ranges and the corresponding
% knot indices along each direction

[elRangeU,elConnU] = buildConnectivity(p,uKnot,noElemsU);
[elRangeV,elConnV] = buildConnectivity(q,vKnot,noElemsV);
[elRangeW,elConnW] = buildConnectivity(r,wKnot,noElemsW);

% combine info from two directions to build the elements
% element is numbered as follows
%  5 | 6 | 7 | 8
% ---------------
%  1 | 2 | 3 | 4
% for a 4x2 mesh

noElems = noElemsU * noElemsV * noElemsW;
element = zeros(noElems,(p+1)*(q+1)*(r+1));

e = 1;
for w=1:noElemsW
    wConn = elConnW(w,:);
    for v=1:noElemsV
        vConn = elConnV(v,:);
        for u=1:noElemsU
            c = 1;
            uConn = elConnU(u,:);
            for i=1:length(wConn)
                for j=1:length(vConn)
                    for k=1:length(uConn)
                        element(e,c) = chan(wConn(i),vConn(j),uConn(k));
                        c = c + 1;
                    end
                end
            end
            e = e + 1;
        end
    end
end

index = zeros(noElems,3);
count = 1;

for i=1:size(elRangeW,1)
    for j=1:size(elRangeV,1)
        for k=1:size(elRangeU,1)
            index(count,1) = k;
            index(count,2) = j;
            index(count,3) = i;
            count = count + 1;
        end
    end
end



jacob   = zeros(3,3);

Xi= zeros(1,p+1);
Eta = zeros(1,q+1);
Zeta = zeros(1,r+1);
dRdxi= zeros(1,p+1);
dRdeta = zeros(1,q+1);
dRdzeta = zeros(1,r+1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 3 ); 

% Assembling system of equation
% Stiffness matrix and external force vector

size_km_elem   = size(element,2);
ASSEMBLY = repmat(struct('Ig', zeros(size_km_elem*size_km_elem,1), ...
                         'Jg', zeros(size_km_elem*size_km_elem,1), ...
                         'Kg', zeros(size_km_elem*size_km_elem,1), ...
                         'Mg', zeros(size_km_elem*size_km_elem,1)), ...
                         noElems,1);

%disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])
%%
% Loop over elements (knot spans)

parfor e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    idw    = index(e,3);
    
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = elRangeW(idw,:); % [zeta_k,zeta_k+1]
 
    sctr   = element(e,:);          %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts sctr+2*noCtrPts]; % scatters a B matrix
    nn     = length(sctr);
    B      = zeros(3,2*nn);
    Ba     = zeros(3,nn);
    Na     = zeros(1,nn);
   
    % loop over Gauss points

    Ks_l = zeros(size_km_elem,size_km_elem);
    Ms_l = zeros(size_km_elem,size_km_elem);
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1)  );
        Eta     = parent2ParametricSpace(etaE,pt(2) );
        Zeta    = parent2ParametricSpace(zetaE,pt(3));
        J2      = jacobianPaPaMapping3d(xiE,etaE,zetaE);
        
        % compute derivative of basis functions w.r.t parameter coord
        
        [N, dRdxi, dRdeta, dRdzeta] = NURBS3DBasisDers([Xi;Eta;Zeta],...
                                   p,q,r,uKnot,vKnot,wKnot,weights');

        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        pts = controlPts(sctr,:);
        
        % Jacobian matrix
             
        jacob = pts'*[dRdxi' dRdeta' dRdzeta'];
        J1    = det(jacob);
        
        % Jacobian inverse and spatial derivatives
        
        %invJacob   = inv(jacob);
        dRdx       = [dRdxi' dRdeta' dRdzeta']/jacob;
        
        Ba(1,1:nn)      = dRdx(:,1)';
        Ba(2,1:nn)      = dRdx(:,2)';
        Ba(3,1:nn)      = dRdx(:,3)';
        Na(1,1:nn)      = N(:)';
     
        % B matrix
        
        B = strainDispMatrix3d(nn,dRdx);
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        % K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * J1 * J2 * wt;
        Ks_l  = Ks_l  + Ba' *  Ba * J1 * J2 * wt;
        Ms_l  = Ms_l  + (1/(c0*c0))*(Na' *  Na) * J1 * J2 * wt;
    end
    ASSEMBLY(e).Ig = kron(ones(size_km_elem,1),sctr');
    ASSEMBLY(e).Jg = kron(sctr',ones(size_km_elem,1));
    ASSEMBLY(e).Kg = reshape(Ks_l,size_km_elem*size_km_elem,1);
    ASSEMBLY(e).Mg = reshape(Ms_l,size_km_elem*size_km_elem,1);
end

Kg = vertcat(ASSEMBLY.Kg);
Mg = vertcat(ASSEMBLY.Mg);
Ig = vertcat(ASSEMBLY.Ig);
Jg = vertcat(ASSEMBLY.Jg);
Ka = sparse(Ig,Jg,Kg,noDofs,noDofs);
Ma = sparse(Ig,Jg,Mg,noDofs,noDofs);
