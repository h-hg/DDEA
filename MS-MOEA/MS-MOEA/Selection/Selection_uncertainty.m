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
function [OffCho,mse]=Selection_uncertainty(Offspring,Offmse,N)
                An = size(Offmse,1);
                Next = zeros(1,N);
                for i = 1:N
                    A_num = randperm(size(Offmse,1));
                    Uncertainty = mean(Offmse(A_num(1:ceil(0.1*An)),:),2);%0.3
                    [~,best]    = max(Uncertainty);
                    Next(i)     = A_num(best);
                end
                OffCho=Offspring(1,Next);
                mse=Offmse(Next,:);

end