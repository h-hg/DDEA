function membership = Paretoset(X)
% PARETOSET  To get the Pareto set from a given set of points.从给定的一组点得到帕累托集。
% synopsis:           membership =paretoset (objectiveMatrix)
% where:
%   objectiveMatrix:目标矩阵 [number of points * number of objectives] array
%   membership:隶属度      [number of points *1] logical vector to indicate if ith
%                    point belongs to the Pareto set (true) or not
%                    (false).用逻辑向量表示该点是否属于pareto集合
%
% by Yi Cao, Cranfield University, 02 June 2007
% Revised by Yi Cao on 17 October 2007
% Version 3, 21 October 2007, new sorting scheme to improve speed.
% Bugfix, 25 July 2008, devided by zero error is fixed.
%
% Examples: see paretoset_examples
%

m=size(X,1);%向量行数
Xmin=min(X);%每列最小值返回为行向量
X1=X-Xmin(ones(m,1),:);     %make sure X1>=0;X中每个元素与Xmin对应元素相减
Xmean=mean(X1); %每列求平均
%sort X1 so that dominated points can be removed quickly对X1进行排序，以便可以快速移除支配点
[~,checklist]=sort(max(X1./(Xmean(ones(m,1),:)+max(Xmean)),[],2));
Y=X(checklist,:);                  
membership=false(m,1);
while numel(checklist)>1
    k=checklist(1);
    [membership(k),checklist,Y]=paretosub(Y,checklist);
end
membership(checklist)=true;

function [ispareto,nondominated,X]=paretosub(X,checklist)

Z=X-X(ones(size(X,1),1),:);
nondominated=any(Z<0,2);                    %retain nondominated points from the check list                         
ispareto=all(any(Z(nondominated,:)>0,2));   %check if current point belongs to pareto set    
X=X(nondominated,:);
nondominated=checklist(nondominated);
