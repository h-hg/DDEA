function [ W,B,C,S] = RBF_EnsembleUN( L,c,nc,T)
% Usage: [ W,B,C,S] = RBF_EnsembleUN( L,c,nc,T)
%Build RBF Model Pool
% Input:
% L             - Offline Data with c Decision Variables and Exact Objective Value
% c             - Number of Decision Variables
% nc            - Number of neurons of RBF models
% T             - %Number of RBF models
%
% Output: 
% W             - Weights of RBF Models
% B             - Bais of RBF Models
% C             - Centers of RBF Models
% S             - Widths of RBF models
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
t=0.5;%probability of out-of-bag
W=zeros(T,nc);
B=zeros(1,T);
C=zeros(c,nc,T);
S=zeros(nc,T);
%Bootstrap Sampling
Index=rand(size(L,1),T);
I1=find(Index<t);
I0=find(Index>=t);
Index(I1)=1;
Index(I0)=0;
%Build Model Pool
for i=1:T
    I= find(Index(:,i)==1);
    L1=L(I,:);
    [ W2,B2,Centers,Spreads ] = RBF( L1(:,1:c),L1(:,c+1),nc);
    W(i,:)=W2;
    B(i)=B2;
    C(:,:,i)=Centers;
    S(:,i)=Spreads;
end

end

