function d = DOEMinDistCriterion(X, varargin)
%% from SRGTSToolbox
%Function DOEMinDistCriterion calculates the minimum distance, dmin,
%criterion for a given design of experiment. The minimum distance criterion
%is a measure of how well spread the points of the sample are over the
%design space:
% 
% dmin = min(D)
%         i
% 
% where D is the matrix of Euclidean distances between two points of the
% design. The general inter-point distance between any point pair can be
% expressed as follows:
% 
%           nv
% d_ij = ( sum |x_ik - x_jk|^(2) )^(1/2)
%          k=1
% 
% where nv is the number of variables.
% 
%Thus, for example:
%
%     dmin = srgtsDOEMinDistCriterion(X): returns the dmin value for the
%     design given in X.
%
%Example:
%     % create Latin hypercube design.
%     NDV     = 2;
%     NPOINTS = 12;
%     LHD = lhsdesign(NPOINTS, NDV)
% 
%     LHD =
% 
%     0.5268    0.8588
%     0.6035    0.3069
%     0.1047    0.5041
%     0.7173    0.1638
%     0.4454    0.2134
%     0.0691    0.7182
%     0.2745    0.4362
%     0.4140    0.9337
%     0.8103    0.4011
%     0.9128    0.0425
%     0.2419    0.6295
%     0.9314    0.7795
% 
%     dmin = srgtsDOEMinDistCriterion(LHD)
% 
%     dmin =
% 
%     0.1354
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
    npoints = length(X(:,1));
else
    npoints = varargin{1};
end

if npoints < 500
    d = pdist(X, 'euclidean');
    d = min(d);
else

    d = Inf;
    for c1 = 1 : (npoints - 1)

        x01 = X(c1,:);

        for c2 = (c1 + 1) : npoints

            x02 = X(c2,:);

            dnew = pdist([x01; x02], 'euclidean');

            if dnew < d
                d = dnew;
            end

        end

    end

end

return
