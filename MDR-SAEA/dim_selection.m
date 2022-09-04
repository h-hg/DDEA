% Input:
% f,g_1,g_2       -Function handle
% sparsity        -Sparsity of the benchmark problem
% D               -Dimension of decision variables
% lower           -Lower bound of decision variables
% upper           -Upper bound of decision variables
%
% Output:
% dim             -Selected out non-zero dimension
% TMask           -New Mask matrix based on selected dimensions
% TempPop         -New population based on selected dimensions
% popobj          -Objective values of TempPop
%%%% Authors: Zheng Tan, Handing Wang, Shulei Liu
%%%% Xidian University, China and Chinese Academy of Military Science, China.
%%%% EMAIL: zhengtan@stu.xidian.edu.cn, hdwang @ xidian.edu.cn
%%%% WEBSITE: https://sites.google.com/site/handingwanghomepage
%%%% DATE: March 2021
% ------------------------------------------------------------------------
% This code is part of the program that produces the results in the following paper:
%
% Zheng Tan, Handing Wang, Shulei Liu, Multi-Stage Dimension Reduction for Expensive Sparse Multi-Objective Optimization Problems, Neurocomputing, vol.440, no.14, pp.159–174, 2021.
%
% You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
% ------------------------------------------------------------------------

function [dim,TMask,TempPop,popobj] = dim_selection(f,sparsity,g_1,g_2,D,lower,upper)

base = zeros(1,D);
Dec = lower+(upper-lower).*rand(D,D);
Mask       = eye(D);
Population = Dec.*Mask;
ind = [base;Population];
popobj_o = f(ind,sparsity,g_1,g_2);
[FrontNo,~] = NDSort(popobj_o,D+1);
dim = find(FrontNo<=FrontNo(1,1))-1;
index = size(dim,2);
t = 1;
TDec_o = Dec;
TMask_o = Mask;
TempPop_o = Population;
popobj_o = popobj_o(2:end,:);
if index<=0.1*D
    while t<3
base = zeros(1,D);
Dec = lower+(upper-lower).*rand(D,D);
Mask       = eye(D);
Population = Dec.*Mask;
ind = [base;Population];
popobj_o = f(ind,sparsity,g_1,g_2);
[FrontNo,~] = NDSort(popobj_o,D+1);
dim_temp = find(FrontNo<=FrontNo(1,1))-1;
dim = union(dim_temp,dim);
t = t + 1;
TDec_o = [TDec_o;Dec];
TMask_o = [TMask_o;Mask];
TempPop_o = [TempPop_o;Population];
popobj_o = [popobj_o;popobj_o(2:end,:)];
    end

else
    while t<3
        base = zeros(1,D);
        Dec = lower+(upper-lower).*rand(D,D);
        Mask       = eye(D);
        Population = Dec.*Mask;
        ind = [base;Population];
        popobj_o = f(ind,sparsity,g_1,g_2);
        indicator = zeros(1,1000);
        for i=1:1000
            indicator(i) = compare(popobj_o(i+1,1),popobj_o(1,1))+compare(popobj_o(i+1,2),popobj_o(1,2));
        end
        index_1 = find(indicator==0);
        index_2 = find(indicator==-1);
        dim_temp = union(index_1,index_2);
        dim = intersect(dim_temp,dim);
        t = t + 1;
        TDec_o = [TDec_o;Dec];
        TMask_o = [TMask_o;Mask];
        TempPop_o = [TempPop_o;Population];
        popobj_o = [popobj_o;popobj_o(2:end,:)];
            end
 end
        dim(dim==0)=[];
        N = size(dim,2);
        TempPop = zeros(3*N,D);
        TempPop(1:N,:) = TempPop_o(dim(1,1:N),:);
        TempPop(N+1:2*N,:) = TempPop_o(D+dim(1,1:N),:);
        TempPop(2*N+1:3*N,:) = TempPop_o(2*D+dim(1,1:N),:);
        TMask = [eye(N);eye(N);eye(N)];
        popobj = zeros(3*N,M);
        popobj(1:N,:) = popobj_o(dim(1,1:N),:);
        popobj(N+1:2*N,:) = popobj_o(D+dim(1,1:N),:);
        popobj(2*N+1:3*N,:) = popobj_o(2*D+dim(1,1:N),:);
end

function i = compare(a,b)
    if a<b
     i = -1;
    end
    if a>b
     i = 1;
    end
    if a==b
     i = 0;
    end
end