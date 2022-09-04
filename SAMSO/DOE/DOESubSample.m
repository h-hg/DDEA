function [SubDOE, SubDOEIdx, RemDOE, RemDOEIdx] = DOESubSample(DOE,npoints,varargin)
%% from SRGTSToolbox
%Function srgtsDOESubSample selects points from a user-defined design
%according to a specific optimality criterion.
% 
% THIS FUNCTION USES THE NATIVE MATLAB FUNCTION: candexch!
% 
%Each row of the design represents one point (or sample). Design variables
%are normalized so that the hypercube points take values between 0 and 1.
%Thus, for example:
%
%     P = srgtsDOESubSample(DOE, NPOINTS): generates an NPOINTS-by-NDV
%     matrix. NPOINTS is the number of points and NDV is the number of
%     variables. It uses the 'MinMax' criterion running for 1000
%     iterations.
%
%     P = srgtsDOESubSample(DOE, NPOINTS, CRITERION, CRITOPT): generates an
%     NPOINTS-by-NDV matrix. NPOINTS is the number of points and NDV is the
%     number of variables. It uses the criterion specified by CRITERION and
%     CRITOPT. In this case:
%          * CRITERION: is a string that defines the method for point
%          selection, and
%          * CRITOPT: is the option used during the point selection.
%
%     CRITERION can be one of the following:
%          * 'MaxMin': select points according to the maximization of the
%          minimum distance between points. In this case, CRITOPT is the
%          number of iterations to be performed.
%          * 'Doptimal': which select points according to the D-optimal
%          criterion, in which designs are obtained by the maximization of
%          the determinant D = |X'X|. This optimality criterion results in
%          minimizing the generalized variance of the parameter estimates
%          for a pre-specified model. As a result, the 'optimality' of a
%          given D-optimal design is model dependent (i.e. it depends on
%          the order of the polynomial response surface). In this case,
%          CRITOPT is a vector in which the first element represents the
%          degree of the polynomial response surface used to model the
%          data; and the second element represents the number of iterations
%          to be performed. (NOT AVAILABLE IN OCTAVE!)
%
%     [P, IDX] = srgtsDOESubSample(...): also returns a vector that
%     informs the indexes of each point of P (i.e., the position of each
%     point of P in the input design, DOE).
%
%     [P, IDX, REMP] = srgtsDOESubSample(...): also returns the remaining
%     points of the input design (i.e., those not used in P).
%
%     [P, IDX, REMP, REMPIDX] = srgtsDOESubSample(...): also returns the
%     vector of indexes of the remaining points of the input design (i.e.,
%     the position of each point of REMP in the input design, DOE).
%
%Example:
%     % create pre-defined design (200x2).
%     NLHPOINTS = 200;
%     NDV       = 2;
%     LHD       = lhsdesign(NLHPOINTS, NDV); % it can be any user-defined design
%
%     % select a design with 16 points usign 'MaxMin'.
%     NPOINTS     = 16;
%     NITERATIONS = 1000;
%     P01         = srgtsDOESubSample(LHD, NPOINTS, 'MaxMin', NITERATIONS)
%
% 
%     P01 =
% 
%     0.6143    0.0262
%     0.7683    0.6416
%     0.5039    0.1999
%     0.8419    0.2663
%     0.3603    0.1220
%     0.6199    0.6631
%     0.8922    0.5801
%     0.0620    0.8628
%     0.0426    0.1108
%     0.4147    0.4425
%     0.0131    0.4033
%     0.2987    0.8845
%     0.6640    0.4400
%     0.0022    0.6471
%     0.8828    0.9007
%     0.5414    0.9683
%     
%     % select a design with 16 points usign 'Doptimal'.
%     PRSDEGREE = 2;
%     P02       = srgtsDOESubSample(LHD, NPOINTS, 'Doptimal', [PRSDEGREE NITERATIONS])
% 
%     P02 =
% 
%     0.9523    0.0728
%     0.0426    0.1108
%     0.9779    0.9616
%     0.9220    0.0503
%     0.0164    0.1728
%     0.0084    0.4736
%     0.9851    0.3844
%     0.4987    0.5545
%     0.2803    0.0046
%     0.0986    0.9899
%     0.9730    0.8878
%     0.9640    0.1583
%     0.0022    0.6471
%     0.4694    0.5324
%     0.1243    0.9745
%     0.4635    0.9812
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
switch nargin
    case 2
        criterion  = 'MaxMin';
        critOption = 1000;
    case 4
        criterion  = varargin{end - 1};
        critOption = varargin{end};
        
    otherwise
        error('SURROGATES Toolbox:srgtsDOESubSample','User must set both "CRITERION" and "CRITOPTIONS," otherwise just enter DOE and NPOINTS to use the default settings of this function.');
end

set = {'MaxMin' 'Doptimal'};
[valid, setstr] = utilitiesIsStringInASet(criterion,set);

if valid == 0
    error('SURROGATES Toolbox:srgtsDOESubSample',sprintf('"CRITERION" must be one of the following: %s.', setstr));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% run point selection

SubDOE    = [];
SubDOEIdx = [];
RemDOE    = [];
RemDOEIdx = [];

[NbTotalPoints NDV] = size(DOE);

switch criterion
    case 'MaxMin'
        [SubDOE, SubDOEIdx, RemDOE, RemDOEIdx] = GetMaxMin(DOE,NbTotalPoints,npoints,NDV,critOption);
    case 'Doptimal'
        [SubDOE, SubDOEIdx, RemDOE, RemDOEIdx] =GetDoptimal(DOE,NbTotalPoints,npoints,NDV,critOption);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friend functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [valid, setstr] = utilitiesIsStringInASet(value,set)
%utilitiesIsStringInASet(VALUE,SET) is true if VALUE is a string in the set.

setstr = sprintf('\n\t* %s',set{:});
valid = 1;

if ~( ischar(value) && any( strcmpi( value , set)) )
    valid = 0;
end    

return
