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
function PopulationObj=cal_objectvalue(PopulationDec,num_obj,choose_Nflag,model)
           if choose_Nflag==1
            for i = 1: size(PopulationDec,1)
                for j = 1 : num_obj
                    [PopulationObj(i,j),~,~] = predictor(PopulationDec(i,:),model{j});
                end
            end
           end
           if choose_Nflag==0
            for j = 1 : num_obj
                PopulationObj(:,j)=cal_via_net(PopulationDec,model{j});
            end
           end
end