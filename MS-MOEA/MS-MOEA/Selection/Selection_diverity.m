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
function Pop=Selection_diverity(population,N)
         
         pop_dec=population.decs;
         pop_obj=population.objs;
         PD_fitness=zeros(1,size(pop_dec,1));
         Location=zeros(1,N);
          total_pd=PD(pop_obj);%目标域
%          total_pd=PD(pop_obj);%决策域
         for i=1:size(pop_dec,1)
%              目标域
           tran_obj=pop_obj;
           tran_obj(i,:)=[];
           PD_fitness(1,i)=total_pd-PD(tran_obj);
%              %决策域
%            tran_dec=pop_dec;
%            tran_dec(i,:)=[];
%            PD_fitness(1,i)=total_pd-PD(tran_dec);
         end
         tran_fitness=PD_fitness;
         for i=1:N

           tran_fitsort=sort(tran_fitness,'descend');
           [~,index]=find(tran_fitness==tran_fitsort(1,i),1);
           Location(i)=index;
           tran_fitness(Location(i))=-1000000;
         end
         for i=1:N
             Pop(1,i)=population(1,Location(i));
         end
end