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
%降噪处理
function [y_deal,T_A,T_A_RBF,T_A_Clu]=data_deal(T_Aori,T_AIori,num_obj,num_vari,choose_Nflag)
        %使用RBF模型进行滤波
          for i = 1 : num_obj
            dmodel     = construct_rbfn(T_Aori(:,1:num_vari),T_Aori(:,num_vari+i),size(T_Aori,1));
            Model{i}   = dmodel;
          end
          for j = 1 : num_obj
            y_RBFdeal(:,j) = cal_via_net(T_Aori(:,1:num_vari), Model{j}); 
          end
           y_RBFdeal=y_RBFdeal;
          %使用聚类思想进行滤波
            y_Cludeal= T_AIori.objs;                                         
          T_A_RBF=[T_Aori(:,1:num_vari),y_RBFdeal];
          T_A_Clu=[T_Aori(:,1:num_vari),y_Cludeal];

            if choose_Nflag==1
                y_deal=y_Cludeal;
                T_A=T_A_Clu;
            else
                y_deal=y_RBFdeal;
                T_A=T_A_RBF;
            end    
                
end                
            

                
                
                