% Input:
% mask_dim        -Mask matrix of obtained population in the evolutionary process
% dim             -Non-zero dimensions obtained in the first stage
% num             -Number of dimensions for the SAEA
%
% Output:
% dim_saea        -Obtained non-zero dimensions in the second stage
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

function dim_saea = dim_base(mask_dim,dim,num)
non_zero = zeros(1,size(mask_dim,1));
for i = 1:size(mask_dim,1)
    non_zero(i) =size(find(mask_dim(i,:)),2);
end
index = find(non_zero==num);%含有num个维度的mask
n = size(index,2);
dim_b = zeros(n,num);
for i = 1:n
    dim_b(i,:) = dim(find(mask_dim(index(i),:)));
end 
% t= tabulate(dim_b);
dim_saea = unique(dim_b)';
end