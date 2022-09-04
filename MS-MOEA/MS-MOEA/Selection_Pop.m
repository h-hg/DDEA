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
function PopDec=Selection_Pop(PopNew_x,PopNew_y,T_AI,Size,kappa)
        T_ADEC=T_AI.decs;
        T_AOBJ=T_AI.objs;
        
        for i=1:size(PopNew_x,1)
            Population(1,i)=INDIVIDUAL(PopNew_x(i,:),PopNew_y(i,:));
        end
        Population(1,size(PopNew_x,1)+1:95) = EnvironmentalSelection(T_AI(1,1:size(T_AI,2)-5),Size-5-size(PopNew_x,1),kappa);   
        PopDec=zeros(Size,size(T_ADEC,2));
        PopDec(1:95,:)=Population(1,1:95).decs;
        non_dominated_set=[];
        while size(non_dominated_set,1)<5
            index = Paretoset(T_AOBJ);   %得到pareto集，找到初始决策点中的pareto前沿点
            set= T_ADEC(index,:);
            T_AOBJ=T_AOBJ(~index,:);
            T_ADEC=T_ADEC(~index,:);
            non_dominated_set=[non_dominated_set;set];
        end
        
        for i=1:5
            for j=1:size(T_ADEC,2)
                non_dominated_set(i,j)= non_dominated_set(i,j)+mvnrnd (0,0.5,1);
                if non_dominated_set(i,j)>1
                    non_dominated_set(i,j)=1;
                elseif non_dominated_set(i,j)<0
                    non_dominated_set(i,j)=0;
                end
            end
        end
         PopDec(96:end,:)=non_dominated_set(1:5,:);

 
end