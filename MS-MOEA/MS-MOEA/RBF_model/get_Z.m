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
function Z = get_Z(popdec, RName,obj_nets)
    kernal_name = RName;
    PopDec = popdec;
    eval(['Z = ', kernal_name, '(PopDec, obj_nets.centers, obj_nets.sigma);']);
    Z = [Z, ones(size(Z, 1), 1)];
end