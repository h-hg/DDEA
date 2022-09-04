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

% This code is written by Nan Zheng.
% Email: nanszheng@163.com
%种群初始化
function [initial_x,initial_y]=Pop_initial(fun_test,Size,num_vari,num_obj,jia)
    PopDec  = lhsamp(Size,num_vari);
    initial_x = repmat(ones(1,num_vari)-zeros(1,num_vari),Size,1).*PopDec+repmat(zeros(1,num_vari),Size,1);
    if jia==1
        if(num_obj==2)
        initial_y = feval(fun_test, initial_x, num_obj)+mvnrnd ([0 0],[0.2 0 ;0 0.2 ],size(initial_x,1));
        else
        initial_y = feval(fun_test, initial_x, num_obj)+mvnrnd ([0 0 0],[0.2 0 0;0 0.2 0;0 0 0.2],size(initial_x,1));
        end
    else
    initial_y = feval(fun_test, initial_x, num_obj);
    end
end