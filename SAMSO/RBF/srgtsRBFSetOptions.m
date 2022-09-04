function srgtOPT = srgtsRBFSetOptions(P, T, FIT_Fn, FIT_LossFn, RBF_type, RBF_c, RBF_usePolyPart)
%Function srgtsRBFSetOptions creates the SURROGATES Toolbox option
%structure for radial basis function models. This structure contains the
%following fiels:
%
%* GENERAL PARAMETERS
%
%   SRGT   - Identifier of the surrogate technique: 'RBF'.
%   P      - NPOINTS-by-NDV matrix, where NPOINTS is the number of points
%            of the sample and NDV is the number of design variables.
%            Default: Empty matrix.
%   T      - NPOINTS-by-1 vector of responses on the P matrix points.
%            Default: Empty vector.
%   FIT_Fn - Function handle of the fitting function (which is used to
%            optimize KRG_Theta). [@rbf_build | @srgtsFit].
%            Default: @rbf_build.
%
%* RADIAL BASIS FUNCTION PARAMETERS
%
% RBF_type        - Type of the basis functions (default = 'MQ'):
%                     'BH' = Biharmonic
%                     'MQ' = Multiquadric
%                     'IMQ' = Inverse Multiquadric
%                     'TPS' = Thin plate spline
%                     'G' = Gaussian
% RBF_c           - Parameter c value (default = 1). Similar to the
%                   "Spread" of native MATLAB radial basis neural network.
% RBF_usePolyPart - Use also the polynomial term P of the model Y = P + RBF
%                   (default = 0, do not use)
%
%The SURROGATES Toolbox uses the RBF toolbox of Jekabsons (2009) to execute
%the radial basis function algorithm.
%
%This is how you can use srgtsRBFSetOptions:
%
%     OPTIONS = srgtsRBFSetOptions: creates a structure with the empty
%     parameters.
%
%     OPTIONS = srgtsRBFSetOptions(P, T): Given the sampled data P (input
%     variables) and T (output variables), it creates a structure with
%     default parameters used for all not specified fields.
%
%     OPTIONS = srgtsRBFSetOptions(P, T, FIT_Fn, RBF_type, RBF_c,
%     RBF_usePolyPart):creates a structure with each of the specified
%     fields.
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
%     options = srgtsRBFSetOptions(X, Y)
%
%     options =
%
%                SRGT: 'RBF'
%                   P: [5x1 double]
%                   T: [5x1 double]
%            RBF_type: 'MQ'
%               RBF_c: 2
%     RBF_usePolyPart: 0
%
%REFERENCES:
%
%Jekabsons G, RBF - Radial Basis Function interpolation for Matlab/Octave,
%Ver. 1.1, Riga Technical University, 2009.
%Available at: http://www.cs.rtu.lv/jekabsons/regression.html.

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
srgtOPT.SRGT = 'RBF';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options
switch nargin
    case 0
        srgtOPT.P = [];
        srgtOPT.T = [];
        
        srgtOPT.FIT_Fn     = [];

        srgtOPT.RBF_type        = [];
        srgtOPT.RBF_c           = [];
        srgtOPT.RBF_usePolyPart = [];
        
    case 2
        srgtOPT.P = P;
        srgtOPT.T = T;
        
        srgtOPT.FIT_Fn     = @rbf_build;

        srgtOPT.RBF_type        = 'MQ';
        srgtOPT.RBF_c           = 2;
        srgtOPT.RBF_usePolyPart = 0;
        
    otherwise
        srgtOPT.P = P;
        srgtOPT.T = T;
        
        srgtOPT.FIT_Fn     = FIT_Fn;

        srgtOPT.RBF_type        = RBF_type;
        srgtOPT.RBF_c           = RBF_c;
        srgtOPT.RBF_usePolyPart = RBF_usePolyPart;
        
end

return
