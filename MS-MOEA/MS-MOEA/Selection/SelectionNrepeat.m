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
function Popcho=SelectionNrepeat(PopNew)
        PopNew1=[];
      for i=1:size(PopNew,1)
          count=0;i=1;
          if size(PopNew,1)==1
             PopNew1(size(PopNew1,1)+1,:)=PopNew;
             break;
          end
          tran_new=PopNew(i,:);
          PopNew(i,:)=[];
          for j=1:size(PopNew,1)
             if norm(tran_new-PopNew(j,:))<=1e-5
                break;
             else
                 count=count+1;
             end             
          end
          if count==size(PopNew,1)
             PopNew1(size(PopNew1,1)+1,:)=tran_new;
          end
      end
      Popcho=PopNew1;

end