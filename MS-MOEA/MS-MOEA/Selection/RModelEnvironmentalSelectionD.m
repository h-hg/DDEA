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
function Pop=RModelEnvironmentalSelectionD(population,MaxSize)
         

         pop_dec=population.decs;
         pop_obj=population.objs;
         p=1/size(pop_obj,2);
         Popobj_Nor_pre = (pop_obj - repmat(min(pop_obj,[],1),size(pop_obj,1),1))...
                        ./repmat(max(pop_obj,[],1) - min(pop_obj,[],1),size(pop_obj,1),1);
         index=Paretoset(pop_obj);
         pop_dec=pop_dec(index,:);
         pop_obj=pop_obj(index,:);
         Popobj_Nor_pre=Popobj_Nor_pre(index,:);
         N=size(pop_obj,1);
         if size(pop_obj,1)<=MaxSize
              Pop=population(1,index);             
             return;
         end
         
            %% Select the extreme solutions first
    
    Choose = false(1,N);
    
    M = size(Popobj_Nor_pre,2);
    select = randperm(M);
     
   

    Choose(select(1)) = true;
    
    
    
    
    

    %% Delete or add solutions to make a total of K solutions be chosen by truncation
    if sum(Choose) > MaxSize
        % Randomly delete several solutions
        Choosed = find(Choose);
        k = randperm(sum(Choose),sum(Choose)-MaxSize);
        Choose(Choosed(k)) = false;
    elseif sum(Choose) < MaxSize
        % Add several solutions by truncation strategy
        Distance = inf(N);
        for i = 1 : N-1
            for j = i+1 : N
                Distance(i,j) = norm(Popobj_Nor_pre(i,:)-Popobj_Nor_pre(j,:),p);
                Distance(j,i) = Distance(i,j);
            end
        end
        while sum(Choose) < MaxSize
            Remain = find(~Choose);
            [~,x]  = max(min(Distance(~Choose,Choose),[],2));
            Choose(Remain(x)) = true;
        end
    end
    Pop=population(1,Choose);
end