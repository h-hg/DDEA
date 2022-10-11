function [MSE, RMSE, RRMSE, R2] = rbftest(model, Xtr, Xtst, Ytst)
% RBFTEST
% Tests an RBF model on a test data set (Xtst, Ytst)
%
% Call
%   [MSE, RMSE, RRMSE, R2] = rbftest(model, Xtr, Xtst, Ytst)
%
% Input
% model     : RBF model
% Xtr       : Inputs of the training data (Xtr(i,:)), i = 1,...,n (the same
%             matrix with which the model was built)
% Xtst, Ytst: Test data points (Xtst(i,:), Ytst(i)), i = 1,...,ntst
%
% Output
% MSE       : Mean Squared Error
% RMSE      : Root Mean Squared Error
% RRMSE     : Relative Root Mean Squared Error
% R2        : Coefficient of Determination

% Gints Jekabsons 2009

if nargin < 4
    error('Too few input arguments.');
end
if (size(Xtst, 1) ~= size(Ytst, 1))
    error('The number of rows in the matrix and the vector should be equal.');
end
if model.n ~= size(Xtr, 1)
    error('The matrix Xtr should be the same matrix with which the model was built.');
end
MSE = mean((rbfpredict(model, Xtr, Xtst) - Ytst) .^ 2);
RMSE = sqrt(MSE);
if size(Ytst, 1) > 1
    RRMSE = RMSE / std(Ytst, 1);
    R2 = 1 - MSE / var(Ytst, 1);
else
    RRMSE = Inf;
    R2 = Inf;
end
return
