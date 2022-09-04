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
function Individual=Evolutionaryperiod(Population,kappa,N)
     
     Next = 1 : length(Population);
     PopulationObj=Population.objs;
    [Fitness,I,C] = CalFitness(PopulationObj,kappa);
    while length(Next) > N-1
        [~,x]   = min(Fitness(Next));
        Fitness = Fitness + exp(-I(Next(x),:)/C(Next(x))/kappa);
        if length(Next)==N
           x1=x;
        end
        Next(x) = [];            
    end
       Individual=Population(1,x1);
end
