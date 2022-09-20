 function [lambda,gamma,c,r]=RFS(sub_size,S,Y,flag,K)
%--------------------------------------------------------------------------
%
%Input:
%sub-size - maximum dimension of all sub-problems  
%S - Sample site matrix (m x d), m=number of sample sites, d=dimension
%Y - Objective function values corresponding to points in S
%flag - string determining what type of RBF model should be used
%K - the number of sub-problems
%
%Output:
%lambda, gamma - vectors with RBF model parameters
%c - the desicion variables of each sub-problem
%r - the samples used for training models of each sub-problem
%--------------------------------------------------------------------------

[n,m]=size(S);

%% RF training
i = 1;
while i<=K
    M = randperm(m,randperm(sub_size,1));  % form sub-problem
    d = length(M);
    L = randperm(n,2*d);
    S1 = S(L,M);     % training set
    Y1 = Y(L,:);
    [lambda{i},gamma{i}] = RBF(S1,Y1,flag);  % train a model for each sub-problem
    c{i} = M;
    r{i} = L;
    if sum(isnan(lambda{i}))==0 
         i=i+1;
    end  
end

end %function