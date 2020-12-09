function [ S ] = SelectModels(W,B,C,S,Xb,c,Q )
% Usage: [ S ] = SelectModels(W,B,C,S,Xb,c,Q )
%Selecting a subset of bagging models
% Input:
% W             - Weights of RBF Models
% B             - Bais of RBF Models
% C             - Centers of RBF Models
% S             - Widths of RBF models
% Xb            - Best Individual
% c             - Number of Decision Variables
% Q             - Number of Selected Models
%
% Output: 
% S             - Index of Selected Models
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

Y = RBF_Ensemble_predictor( W,B,C,S,Xb(:,1:c),c );
T=size(Y,2);
[A,I]=sort(Y);

S=[1:Q]';
S=ceil(S*T/Q);
S=I(S)';

end

