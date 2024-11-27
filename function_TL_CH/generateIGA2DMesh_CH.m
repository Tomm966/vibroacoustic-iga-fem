function [element,elRangeU,elConnU,elRangeV,elConnV,noElems,index] = generateIGA2DMesh_CH(uKnot, ...
                                          vKnot, ...
                                          noPtsX,...
                                          noPtsY, ...
                                          p,q)
% WRITTEN BY TOMMASO LANDI AND CHRISTOPHE HOAREAU

% INPUT:
% uKnot: Knot vector for the xi direction.
% vKnot: Knot vector for the eta direction.
% noPtsX: Number of control points in the xi direction.
% noPtsY: Number of control points in the eta direction.
% p: Degree of the NURBS in the xi direction.
% q: Degree of the NURBS in the eta direction.
% OUTPUT:
% element: A matrix detailing the connectivity of the elements.
% elRangeU: Ranges of the elements in the xi direction.
% elConnU: Connectivity in the xi direction.
% elRangeV: Ranges of the elements in the eta direction.
% elConnV: Connectivity in the eta direction.
% noElems: Total number of elements.
% index: An index matrix representing the location of each element in the mesh.

    uniqueUKnots   = unique(uKnot);
    uniqueVKnots   = unique(vKnot);
    
    noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
    noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % chan = 
    %        1 2 3 4
    %        5 6 7 8
    %        9 10 11 12
    % for a 4x3 control points
    
    
    chan           = zeros(noPtsY,noPtsX);
    
    count = 1;
    
    for i=1:noPtsY
        for j=1:noPtsX
            chan(i,j) = count;
            count = count + 1;
        end
    end
    
    % determine our element ranges and the corresponding 
    % knot indices along each direction
    
    [elRangeU,elConnU] = buildConnectivity(p,uKnot,noElemsU);
    [elRangeV,elConnV] = buildConnectivity(q,vKnot,noElemsV);
    
    % combine info from two directions to build the elements
    % element is numbered as follows
    %  5 | 6 | 7 | 8
    % ---------------
    %  1 | 2 | 3 | 4 
    % for a 4x2 mesh
    
    noElems = noElemsU * noElemsV;
    
    element = zeros(noElems,(p+1)*(q+1));
    
    e = 1;
    for v=1:noElemsV
        vConn = elConnV(v,:);
        for u=1:noElemsU
            c = 1;
            uConn = elConnU(u,:);
            for i=1:length(vConn)
                for j=1:length(uConn)
                  element(e,c) = chan(vConn(i),uConn(j));
                  c = c + 1;
                end
            end
            e = e + 1;
        end
    end
    
    index = zeros(noElems,2);
    count = 1;
    
    for j=1:size(elRangeV,1)
        for i=1:size(elRangeU,1)
            index(count,1) = i;
            index(count,2) = j;   
            count = count + 1;
        end
    end

end