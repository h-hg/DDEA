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
function Pop=ModelSelection_uncertainty(population,Archive,N)
        
             Popdec=population.decs;
             Adec=Archive.decs;
             
            % º∆À„æ‡¿Î
         for i=1:size(Popdec,1)
             for j=1:size(Adec,1)
                 dis(i,j)=norm(Popdec(i,:)-Adec(j,:));
             end
         end
          %% æ‡¿Î≈≈–Ú
         for i=1:size(Popdec,1)
             dis(i,:)=sort(dis(i,:)); 
         end
         tran=dis(:,1);
         tran_sort=sort(dis(:,1),'descend');
         Location=zeros(1,3);
         for i=1:N
             Location(i)=find(tran==tran_sort(i),1);
             tran(Location(i))=0;
             Pop(1,i)=population(1,Location(i));
         end
end