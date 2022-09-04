% Input:
% f,g_1,g_2       -Function handle
% dim_saea        -Obtained non-zero dimensions
% sparsity        -Sparsity of the benchmark problem
% num_temp        -Obtained number of non-zero dimensions
% D               -Dimension of decision variables
% lower           -Lower bound of decision variables
% upper           -Upper bound of decision variables
%
% Output:
% dim             -Selected out non-zero dimensions for the SAEA
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
function dim = dim_final_selection(dim_saea,f,sparsity,g_1,g_2,num_temp,D,lower,upper)
    N = size(dim_saea,2);
    Dec = lower+(upper-lower).*rand(N,D);
        Mask = zeros(N,D);
    for i =1:N
        Mask(i,dim_saea(i)) = 1;
    end
    Population = Dec.*Mask;
    popobj = f(Population,sparsity,g_1,g_2);
    [FrontNo,~] = NDSort(popobj,N);
    [~,I] = sort(FrontNo,'ascend');
    dim = dim_saea(I(1:num_temp));
end