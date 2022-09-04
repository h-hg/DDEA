function phip =DOEPHIpCriterion(X, p, t)
%% from SRGTSToolbox
%Function srgtsDOEPHIpCriterion calculates the PHIp criterion for a given
%design of experiment. The PHIp criterion is a measure of how well spread
%the points of the sample are over the design space:
% 
%           s
% PHIp = ( sum J_i d_i^(-p) )^(1/p)
%          i=1
% 
% where p is a positive integer  d_i's are distance values; J_i is the
% number of point pairs in the design separated by d_i; and s is the number
% of distinct distance values. The general inter-point distance between any
% point pair can be expressed as follows:
% 
%           nv
% d_ij = ( sum |x_ik - x_jk|^(t) )^(1/t)
%          k=1
% 
% where nv is the number of variables.
% 
%Thus, for example:
%
%     PHIP = srgtsDOEPHIpCriterion(X): returns the PHIp value for the
%     design given in X; assuming p = 50 and t = 1.
%
%     PHIP = srgtsDOEPHIpCriterion(X, p): returns the PHIp value for the
%     design given in X with p and assuming t = 1.
%
%     PHIP = srgtsDOEPHIpCriterion(X, p, t): returns the PHIp value for the
%     design give in X with p and t given values.
%
%Example:
%     % create Latin hypercube design.
%     NDV     = 2;
%     NPOINTS = 12;
%     LHD = lhsdesign(NPOINTS, NDV)
% 
%     LHD =
% 
%     0.1284    0.6292
%     0.8696    0.1489
%     0.4507    0.2768
%     0.9840    0.3851
%     0.6049    0.8745
%     0.7308    0.7117
%     0.2591    0.1884
%     0.2403    0.4996
%     0.0263    0.9979
%     0.3747    0.7539
%     0.7752    0.5579
%     0.5071    0.0718
% 
%     phip = srgtsDOEPHIpCriterion(LHD)
% 
%     phip =
% 
%     5.0435
%     
%Results may change from run to run because of the random nature of the
%Latin hypercube design.

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
% check inputs
if nargin < 2
    p = 50;
end
if nargin < 3
    t = 1;
end

[np ndv] = size(X);

phip = 0;
for c1 = 1 : np

    x01 = X(c1,:);

    for c2 = (c1 + 1) : np

        x02  = X(c2,:);
        d    = norm((x01 - x02), t);
        phip = phip + d^(-p);

    end

end

phip = phip^(1.0/p);

return
