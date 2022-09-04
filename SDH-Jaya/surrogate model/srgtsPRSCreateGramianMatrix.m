function X = srgtsPRSCreateGramianMatrix(x, NbVariables, PRSDegree, PRSRemovedIdx)
%Function srgtsPRSCreateGramianMatrix generates the Gramian matrix. Thus, for
%example:
%
%     GRAMIAN = srgtsPRSCreateGramianMatrix(X, NDV, DEGREE): return
%     respective the GRAMIAN of the matrix X, according to the number of
%     design variables given by NDV and the DEGREE of the polynomial
%     response surface.
%
%     GRAMIAN = srgtsPRSCreateGramianMatrix(X, NDV, DEGREE, PRSREMOVEDIDX):
%     performs the same calculations however it eliminates the columns
%     specified in PRSREMOVEDIDX.
%
%Example:
%     % basic information about the problem
%     myFN = @cos;  % this could be any user-defined function
%     designspace = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     % create DOE
%     npoints    = 5;
%     nvariables = 1;
%     X = linspace(designspace(1), designspace(2), npoints)';
%
%     % create the Gramian of X
%     PRSdegree  = 2;
%     GRAMIAN = srgtsPRSCreateGramianMatrix(X, nvariables, PRSdegree)
% 
%     GRAMIAN =

%     1.0000         0         0
%     1.0000    1.5708    2.4674
%     1.0000    3.1416    9.8696
%     1.0000    4.7124   22.2066
%     1.0000    6.2832   39.4784

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NbPoints  = length(x(:,1));

% create X
X(:,1) = ones(NbPoints,1);

if PRSDegree > 0
    X      = [X x];
    X1     = x;
    nc1    = NbVariables;
    n_loc  = [1:NbVariables];
    n_loc1 = 1;
    for i = 2 : PRSDegree
        [nr,nc] = size(X1);
        X2      = [];
        ctr     = 1;
        for k = 1 : NbVariables
            l_ctr = 0;
            for j = n_loc(k) : nc

                X2(:,ctr) = x(:,k).*X1(:,j);
                ctr       = ctr+1;
                l_ctr     = l_ctr+1;
            end
            n_loc1(k+1) = l_ctr + n_loc1(k);
        end
        nc1   = nc;
        X     = [X X2];
        X1    = X2;
        n_loc = n_loc1;
    end

end

if nargin == 4
    X( : , PRSRemovedIdx ) = [];
end

return
