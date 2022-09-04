function [SubDoE, SubDoEIdx, RemDoE, RemDoEIdx] = GetDoptimal(DoE,npointstotal,npoints, ndv, critOption)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PRSdegree   = critOption(1);
NIterations = critOption(2);

SubDoE    = [];
SubDoEIdx = [];
RemDoE    = [];
RemDoEIdx = [];

xmin = zeros(ndv,1);
xmax = ones(ndv,1);

xy_org = DoE;
xy     = DoE;
[X]    = PRSCreateGramianMatrix(xy, ndv, PRSdegree);

r = candexch(X, npoints, 'maxiter', NIterations, 'display', 'off');

% Check for duplicate points
% Array 1 (size npointstotal) has numbers 0 or 1: 0 means point is not used 1
% means point is used
lhs_pts = zeros(npointstotal,1);
r_arr = zeros(npoints,1);
% Array 2 (size npoints) has number 0 (no duplicates) or 1 (duplicates)
for j=1:npoints
    if(lhs_pts(r(j))==0)
        lhs_pts(r(j))=1;
    else
        r_arr(j) = 1;
    end
end

% Remove duplicates from neighboring points (use +/- 5 points un-used
% points in second pass to swap with duplicates which is closest
dup_arr = find(r_arr > 0);
for j=1:length(dup_arr)
    % Construct X matrix for computing determinant
    X_a = X(r(:),:);
    lhs_free = find(lhs_pts == 0);
    % Pick 5 points which are not used and close by
    d_pt = r(dup_arr(j));

    % Compute the distance between this point and all available options
    clear cand_arr
    for k=1:length(lhs_free)
        X_t = X_a;
        cpt = lhs_free(k);
        X_t(j,:)=X(cpt,:);
        cand_arr(k,1)= cpt;
        cand_arr(k,2) = det(X_t'*X_t);
    end

    % Pick minimum distance point
    [mn_v,idx] = max(cand_arr(:,2));
    r(dup_arr(j)) = cand_arr(idx,1);
    lhs_pts(cand_arr(idx,1))=1;
end

SubDoEIdx(:,1) = r(1:npoints);
SubDoE     = xy( SubDoEIdx , : );

RemDoEIdx = [1:npointstotal]';
for counter = 1 : npoints
    idx = find( RemDoEIdx == SubDoEIdx(counter) );
    RemDoEIdx = RemDoEIdx([ [1:idx-1] , [(idx + 1) : length(RemDoEIdx) ] ]);
    
end

RemDoE = DoE( RemDoEIdx , : );

return
