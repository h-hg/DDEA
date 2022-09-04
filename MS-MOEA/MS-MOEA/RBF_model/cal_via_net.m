% An Adaptive Model Switch-based Surrogate-Assisted Evolutionary Algorithm 
% for Noisy Expensive Multi-Objective Optimization
%------------------------------- Reference --------------------------------
% N. Zheng, H. Wang, and B. Yuan, An Adaptive Model Switch-based 
% Surrogate-Assisted Evolutionary Algorithm for Noisy Expensive 
% Multi-Objective Optimization, Complex & Intelligent Systems, accepted, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 HandingWangXD Group. Permission is granted to copy and
% use this code for research, noncommercial purposes, provided this
% copyright notice is retained and the origin of the code is cited. The
% code is provided "as is" and without any warranties, express or implied.

% This code is written by Zhening Liu.
function Offspring = cal_via_net(Decs, obj_nets)
    %%
    [N, ~] = size(Decs);
    Objs = zeros(N, 1);
    %% calculate for evr model && return the average predict
    kernal_name = obj_nets.name;
    eval(['Z = ', kernal_name, '(Decs, obj_nets.centers, obj_nets.sigma);']);
    Z    = [Z, ones(size(Z, 1), 1)];
    Objs = Z*obj_nets.weight+Objs;
    %%
    Offspring = Objs;
end