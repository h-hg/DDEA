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
function h = cal_Convergence(PopObj1,PopObj2,Zmin,p)
% Select part of the solutions in the last front
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    if N1~=N2

       h = 0;
    else
    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    PopObj = (PopObj)./repmat(max(PopObj,[],1) - Zmin,size(PopObj,1),1);
%     [N,M]  = size(PopObj);

    Distance1 = zeros(1,N1);
    Distance2 = zeros(1,N2);
    
%       Extreme = zeros(1,M);
%     w       = zeros(M)+1e-6+eye(M);
%     for i = 1 : M
%         [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
%     end
%     Hyperplane = PopObj(Extreme,:)\ones(M,1);
%     a = 1./Hyperplane;
%     if any(isnan(a))
%         a = max(PopObj,[],1)';
%     end
%     % Normalization
%     PopObj = PopObj./repmat(a',N,1);
    
    
    for i = 1:N1
%     Distance1(i) = norm(PopObj(i,:),p);
    Distance1(i) = sqrt(sum(PopObj(i,:),2));
    end
    
    
    
    for i = 1: N2
%       Distance2(i) = norm(PopObj(N1+i,:),p);
    Distance2(i) = sqrt(sum(PopObj(N1+i,:),2));
    end
    [~,h,~,r1,r2]=signrank_new(Distance1, Distance2,'alpha',0.05);%Wilcoxon符号秩检验，检验两组数据的相关性无需知道数据的分布情况
%     [~,h]  =  ranksum( Distance1, Distance2,'alpha',0.05);
    if h == 1 && (r1-r2 <0)
            h = 0;
    end
%      [p,h,stats,r1,r2]=signrank_new(Distance1, Distance2);
    end
end