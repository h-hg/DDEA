function srgtOPT = srgtsPRSSetOptions(P, T, PRS_Degree, PRS_Regression)
%Function srgtsPRSSetOptions creates the SURROGATES Toolbox option
%structure for polynomial response surface models. This structure contains
%the following fields:
%
%* GENERAL PARAMETERS
%
%   SRGT - Identifier of the surrogate technique: 'PRS'.
%   P    - NPOINTS-by-NDV matrix, where NPOINTS is the number of points of
%          the sample and NDV is the number of design variables.
%          Default: Empty matrix.
%   T    - NPOINTS-by-1 vector of responses on the P matrix points.
%          Default: Empty vector.
%
%* POLYNOMIAL RESPONSE SURFACE PARAMETERS
%
%   PRS_Degree     - Degree of the full polynomial surface that fits
%                    the data. [ positive integer ]. Default: 2.
%   PRS_Regression - Regression technique used to build the PRS
%                    model. [ string | 'Full' | 'StepwiseSRGTS' |
%                    'StepwiseMATLAB' | 'ZeroIntercept' ].
%                    Default: 'Full'.
%                        * Full           : full polynomial (all terms)
%                        * StepwiseSRGTS  : uses an in-house stepwise code.
%                        * StepwiseMATLAB : uses the native MATLAB function
%                                           @stepwisefit with all
%                                           coefficients initially
%                                           "inmodel".
%                        * ZeroIntercept  : similar to 'Full' but with the
%                                           constant term set to zero.
%
%This is how you can use srgtsPRSSetOptions:
% 
%     OPTIONS = srgtsPRSSetOptions: creates a structure with the empty
%     parameters.
%
%     OPTIONS = srgtsPRSSetOptions(P, T): Given the sampled data P (input
%     variables) and T (output variables), it creates a structure with
%     default parameters used for all not specified fields.
%
%     OPTIONS = srgtsPRSSetOptions(P, T, PRS_Degree, PRS_Regression): 
%     it creates a structure with each of the specified fields.
%
%Example:
%     % basic information about the problem
%     myFN = @cos;  % this could be any user-defined function
%     designspace = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     % create DOE
%     npoints = 5;
%     X = linspace(designspace(1), designspace(2), npoints)';
%
%     % evaluate analysis function at X points
%     Y = feval(myFN, X);
%
%     % fit surrogate models
%     options = srgtsPRSSetOptions(X, Y)
%
%     options =
%
%               SRGT: 'PRS'
%                  P: [5x1 double]
%                  T: [5x1 double]
%         PRS_Degree: 2
%     PRS_Regression: 'Full'

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
% data
srgtOPT.SRGT = 'PRS';

srgtOPT.P = P;
srgtOPT.T = T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options
if nargin == 2
    srgtOPT.PRS_Degree     = 2;
    srgtOPT.PRS_Regression = 'Full';
else
    srgtOPT.PRS_Degree     = PRS_Degree;
    srgtOPT.PRS_Regression = PRS_Regression;
end

return
