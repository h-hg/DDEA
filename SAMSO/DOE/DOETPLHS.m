function LHS = DOETPLHS(nPoints, nDV, seed)
%% from SRGTSToolbox
%Function srgtsDOETPLHS generates Latin hypercube design by using the
%Translational Propagation algorithm (TPA). The objective is obtaining
%optimal (or near optimal) Latin hypercube designs without using formal
%optimization. This procedure requires minimal computational effort with
%results virtually provided in real time. The algorithm exploits patterns
%of point locations for optimal Latin hypercube designs based on the
%PHIp-criterion (a variant of the maximum distance criterion). Small
%building blocks (called SEEDs), consisting of one or more points each, are
%used to recreate these patterns by simple translation in the hyperspace.
%Research conducted during development of the TPA found that (i) the
%distribution of PHIp tends to lower values as the dimensionality is
%increased; and (ii) the Latin hypercube designs, obtained via TPA,
%represent an attractive alternative to the optimum Latin hypercube designs
%up to medium dimensions. It is concluded that for up to six dimensions (no
%matter the point density) the proposed Latin hypercube designs provide
%computationally inexpensive estimates of the optimal Latin hypercube
%design.
%Each row of the design represents one run (or sample). Design variables
%are normalized so that the hypercube points take values between 0 and 1.
%Thus, for example:
%
%     P = srgtsDOETPLHS(NPOINTS, NDV): generates an NPOINTS-by-NDV
%     matrix. NPOINTS is the number of points and NDV is the number of
%     variables. The seed design used in this case is a single point placed
%     at the origin of the design space.
%
%     P = srgtsDOETPLHS(NPOINTS, NDV, SEED): generates an NPOINTS-by-NDV
%     matrix. NPOINTS is the number of points and NDV is the number of
%     variables. SEED is the basic Latin hypercube design used to build the
%     ELHD. No normalization is needed for a 1-by-NDV SEED.
%
%     P = srgtsDOETPLHS(NPOINTS, NDV, NTRIALS): generates an NPOINTS-by-NDV
%     matrix. NPOINTS is the number of points and NDV is the number of
%     variables. The algorithm is run NTRIALS times, with seed size varying
%     from 1 to NTRIALS. P is the best design found in terms of the PHIp
%     criterion. The PHIp criterion is a measure of how well spread the
%     points of the sample are over the design space. See
%     srgtsDOEPHIpCriterion for more details about the PHIp calculation.
%
%Example:
%     % create a 14x2 design.
%     NPOINTS = 14;
%     NDV     = 2;
%
%     P = srgtsDOETPLHS(NPOINTS, NDV)
% 
%     P =
% 
%     0.0000    0.0000
%     0.2667    0.0667
%     0.5333    0.1333
%     0.8000    0.2000
%     0.0667    0.2667
%     0.3333    0.3333
%     0.6000    0.4000
%     0.8667    0.4667
%     0.1333    0.5333
%     0.4000    0.6000
%     0.6667    0.6667
%     0.9333    0.7333
%     0.2000    0.8000
%     0.4667    0.8667
%     0.7333    0.9333
%     1.0000    1.0000
%
%For further reference, see:
%Viana FAC, Venter G, and Balabanov V, "An algorithm for fast optimal Latin
%hypercube design of experiments," International Journal for Numerical
%Methods in Engineering, Vol. 82 (2), pp. 135-156, 2010 (DOI:
%10.1002/nme.2750).

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
if nargin == 2
    seed = ones(1, nDV);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
LHS = zeros(nPoints, nDV); % just for memory allocation

if utilitiesIsScalar(seed)
    LHS    = cell(1, seed);
    LHSbig = cell(1, seed);
    phip   = zeros(1, seed);
    for c1 = 1 : seed
        newseed = buildTPLHS(c1, nDV, ones(1,nDV));
        [LHSaux LHSbigaux] = buildTPLHS(nPoints, nDV, newseed);
        
        LHS{c1}    = LHSaux;
        LHSbig{c1} = LHSbigaux;
        phip(c1)   =DOEPHIpCriterion(LHSaux);
    end
    
    idx = find(phip == min(phip)); 
    
    LHS    = LHS{idx(1)};
    
else
    LHS = buildTPLHS(nPoints, nDV, seed);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friend functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = utilitiesIsScalar(s)

[nbr nbc] = size(s);

out = ( (nbr == 1) && (nbc == 1) ) && isa(s,'double');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% build TPLHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [LHS LHSbig] = buildTPLHS(nPoints, nDV, seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% determine the number of points and size of seed and block
nPointsSeed = length( seed(:,1) );
nBlocks     = nPoints/nPointsSeed;
nDiv        = nBlocks^( 1/nDV );

nBlocksSLHBig = nBlocks;
nDivSLHBig    = ceil( nDiv );        % this will create the bigger LHS
if ( nDiv ~= nDivSLHBig );           % need to create a bigger LHS
    nBlocksSLHBig = nDivSLHBig^nDV;
end

nPointsSLHBig = nBlocksSLHBig*nPointsSeed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% reshape seed in order to fit into the current LHBig

if nPointsSeed == 1
    seed = ones(1, nDV); % arbitrarily put at the origin
else
    seed = reshapeSeed( seed, nDV, nPointsSLHBig, nDivSLHBig , nPointsSeed );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% create LHSbig
LHSbig = seed;
displacement = ones(1, nDV); % just for allocation
for c1 = 1 : nDV % move towards one direction at a time

    seed = LHSbig; % update seed with the latest points added

    displacement(1 : (c1 - 1))   = nDivSLHBig^(c1 - 2);
    displacement(c1)             = nPointsSLHBig/nDivSLHBig;
    displacement((c1 + 1) : end) = nDivSLHBig^(c1 - 1);

    for c2 = 2 : nDivSLHBig % fill each of the divisions

        seed   = translateSeed(seed, displacement, nDV);
        LHSbig = vertcat(LHSbig, seed);

    end

end

% shrink ELH if necessary
if (nPoints ~= nPointsSLHBig)
    LHS = shrinkSLH(LHSbig, nPoints, nPointsSLHBig, nDV);
else
    LHS = LHSbig;
end

% scale LHS and LHSbig within the interval [0 1]
NormalizedSpace = [zeros(1,nDV); ones(1,nDV)];
SLHSpace        = [ones(1,nDV); nPoints*ones(1,nDV)];
SLHBigSpace     = [ones(1,nDV); nPointsSLHBig*ones(1,nDV)];

LHS    = srgtsScaleVariable(LHS,    SLHSpace,    NormalizedSpace);
LHSbig = srgtsScaleVariable(LHSbig, SLHBigSpace, NormalizedSpace);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% reshape seed to the current size of the DOE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function seed = reshapeSeed( seed, nDV, nPointsSLHBig, nDivSLHBig , nPointsSeed )

b = ( (nPointsSLHBig / nDivSLHBig) - nDivSLHBig*(nDV - 1) + 1 );

seed = srgtsScaleVariable(seed, [zeros(1,nDV); ones(1,nDV)], [ones(1,nDV); b*ones(1,nDV)]);
seed = round(seed); % just to make sure that the numbers are integer

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% move seed (in each step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function seed = translateSeed(seed, displacement, nDV)

nPoints = length(seed(:,1));

for c1 = 1 : nPoints
    seed(c1,:) = seed(c1,:) + displacement;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% shrink a big TPLHS to the size that was specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function LHS = shrinkSLH(LHSbig, nPoints, nPointsSLHBig, nDV)

% center   = ones(1,nDV)/2;
center   = nPointsSLHBig*ones(1,nDV)/2;
distance = getDistance(LHSbig, center, nPointsSLHBig);

[dummy, idx] = sort(distance);
LHS = LHSbig( idx(1:nPoints), : ); % shrink to nPoints

% now we need to re-establish the LH conditions
minSLH = min(LHS);
for c1 = 1 : nDV
    % place LHS in the origin
    LHS = sortrows(LHS, c1);
    LHS(:,c1) = LHS(:,c1) - minSLH(c1) + 1;

    % eliminate empty coordinates
    flag = 0;
    while flag == 0;
        mask = (LHS(:,c1) ~= ([1:nPoints]'));
        flag = isequal(mask,zeros(nPoints,1));
        LHS(:,c1) = LHS(:,c1) - (LHS(:,c1) ~= ([1:nPoints]'));
    end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% calculate the distance between points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function distance = getDistance(LHS, center, nPoints)

distance = zeros(nPoints,1); % memory allocation

for c1 = 1 : nPoints
    distance(c1) = norm( ( LHS(c1,:) - center) );
end

return
