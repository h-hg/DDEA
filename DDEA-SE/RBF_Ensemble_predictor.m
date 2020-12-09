function [ Y ] = RBF_Ensemble_predictor( W,B,C,S,U,c )
% Usage: [ Y ] = RBF_Ensemble_predictor( W,B,C,S,U,c )
%RBF Predictors 
% Input:
% W             - Weights of RBF Models
% B             - Bais of RBF Models
% C             - Centers of RBF Models
% S             - Widths of RBF models
% U             - Test Data with c Decision Variables
% c             - Number of Decision Variables
%
% Output: 
% Y             - Predictions of RBF Models for U
%
    %%%%    Authors:    Handing Wang, Yaochu Jin, Chaoli Sun, John Doherty
    %%%%    University of Surrey, UK and Taiyuan University of Science and Technology, China.
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage
    %%%%    DATE:       May 2018
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:

%Handing Wang, Yaochu Jin, Chaoli Sun, John Doherty, Offline data-driven evolutionary optimization using selective surrogate ensembles, IEEE Transactions on Evolutionary Computation, Accepted.

%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
Y=[];
T=size(B,2);%Number of RBF models
for i=1:T
    Yhat=RBF_predictor(W(i,:),B(i),C(:,:,i),S(:,i),U(:,1:c));
    Y=[Y,Yhat];
end

end

