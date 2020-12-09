function [TestNNOut]=RBF_predictor(W2,B2,Centers,Spreads,TestSamIn)
% Usage: [TestNNOut]=RBF_predictor(W2,B2,Centers,Spreads,TestSamIn)
% Single RBF Predictor
% Input:
% W2             - Weights of RBF Model
% B2             - Bais of RBF Model
% Centers        - Centers of RBF Model
% Spreads        - Widths of RBF model
% TestSamIn      - Test Data

%
% Output: 
% TestNNOut      - Prediction of RBF Model for TestSamIn
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
N=size(TestSamIn,1);
TestDistance = dist(Centers',TestSamIn');
TestSpreadsMat = repmat(Spreads,1,N);
TestHiddenUnitOut = radbas(TestDistance./TestSpreadsMat);
TestNNOut = W2*TestHiddenUnitOut+B2;
TestNNOut=TestNNOut';
end