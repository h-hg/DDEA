function [yhat predvar] = srgtsPRSPredictor(x, P, srgtSRGT)
%Function srgtsPRSPredictor returns the predicted response and the
%estimated prediction variance of a polynomial response surface model.
%Thus, for example:
%
%     [YHAT PREDVAR] = srgtsPRSPredictor(X, P, SURROGATE): returns both the
%     predicted response YHAT and the prediction variance PREDVAR of the
%     polynomial response surface model SURROGATE at all X sites. P is the
%     experimental design matrix. X can be either a single row vector
%     (single point, with each column representing a variable) or a matrix
%     (each row represents a point). YHAT and PREDVAR are NPOINTS-by-1
%     vectors.
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
%     options = srgtsPRSSetOptions(X, Y);
% 
%     surrogate = srgtsPRSFit(options);
%
%     % create test points
%     Xtest = linspace(designspace(1), designspace(2), 100)';
%
%     % evaluate surrogate at Xtest
%     [Yhat PredVar] = srgtsPRSPredictor(Xtest, X, surrogate);
%
%     plot(X, Y, 'o', ...
%          Xtest, Yhat, ...
%          Xtest, Yhat + sqrt(PredVar), 'r', ...
%          Xtest, Yhat - sqrt(PredVar), 'r')

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

% Gramian matrices
X    = srgtsPRSCreateGramianMatrix(x, srgtSRGT.NbVariables, srgtSRGT.PRS_Degree, srgtSRGT.PRS_RemovedIdx);
X0   = srgtsPRSCreateGramianMatrix(P, srgtSRGT.NbVariables, srgtSRGT.PRS_Degree, srgtSRGT.PRS_RemovedIdx);
Xinv = pinv(X0'*X0);

yhat = X*srgtSRGT.PRS_Beta;

% square of the standard error
stderr2 = srgtSRGT.PRS_SE^2;

NbPointsTest = length(x(:,1));
predvar      = zeros(NbPointsTest, 1);
for c1 = 1 : NbPointsTest
    x = X(c1,:)'; % for proper representation
    predvar(c1,1) = stderr2*(x')*Xinv*x;
end

return
