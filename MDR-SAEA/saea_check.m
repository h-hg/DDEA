% Input:
% dim             -Obtained non-zero dimensions 
% num             -Obtained number of non-zero dimensions for the SAEA
%
% Output:
% result          -Flag of whether the traversal search can be performed
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

% coding:utf-8
function result = saea_check(dim,num)
    global evaluation;
    temp = size(dim,2);
    if 40*nchoosek(temp,num)*num>evaluation
        result = 1;
    else
        result = -1;
    end
end