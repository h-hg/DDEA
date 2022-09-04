function [srgtSRGT srgtSTT] = srgtsRBFFit(srgtOPT)
%Function srgtsRBFFit fits the specified radial basis function model using
%the RBF toolbox of Jekabsons (2009).
%
%    [srgtSRGT srgtSTT] = srgtsRBFFit(srgtOPT)
%
%srgtSRGT is the surrogate structure that contains the following fields:
%* RBF_Model: RBF model, a struct with the elements:
%    * n      : Number of data points in the training data set
%    * meanY  : Mean of Ytr
%    * bf_type: Type of the basis functions
%    * bf_c   : Parameter c value
%    * poly   : Use also the polynomial term
%    * coefs  : Coefficients of the model
%
%srgtSTT is the state structure that contains the following fields:
%* FIT_Fn       : function handle of the fitting function.
%* FIT_FnVal    : value of the loss function (after fitting).
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
%     options = srgtsRBFSetOptions(X, Y);
%
%     surrogate = srgtsRBFFit(options)
%
%     surrogate =
%
%             P: [5x1 double]
%     RBF_Model: [1x1 struct]
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

switch func2str(srgtOPT.FIT_Fn)
    case 'rbf_build'
        srgtSRGT.P = srgtOPT.P;
        srgtSTT = srgtsFitCreateState(srgtOPT);
        srgtSRGT.RBF_Model = rbf_build(srgtOPT.P, srgtOPT.T, ...
            srgtOPT.RBF_type, srgtOPT.RBF_c, srgtOPT.RBF_usePolyPart, 0);
        
    case 'srgtsXVFit'
        [srgtSRGT srgtSTT] = srgtsXVFit(srgtOPT);
        
end

return
