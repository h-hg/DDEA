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

% This code is written by Zhening Liu.
function centers = get_center(popdec,center_num)
%%%Get center points via Kmeans
    %%
    PopDec  = popdec;
    [N, D]  = size(PopDec);
    [~,centers] = kmeans(PopDec,center_num);
%     %LHS to get initial centers
%     centers = zeros(center_num, D);
%     for j=1:D
%         centers(:, j) = linspace(Global.lower(j), Global.upper(j), center_num)';
%     end
%     %´òÂÒ¸÷Î¬¶ÈË³Ðò
%     for j=1:D
%         ran_index     = randperm(size(centers, 1));
%         centers(:, j) = centers(ran_index, j);
%     end
%     
%     %%
%     %Kmeans start
%     index_box   = zeros(N, 1);
%     new_centers = centers;
%     rol_num = 1;
%     while rol_num<=5000
%         [~, index_box]   = min(pdist2(PopDec, centers),[],2);
%         
%         for i=1:center_num
%             center_index      = find(index_box == i);
%             if isempty(center_index)
%                 randnum = randi([1, D*Global.san_den-1], 1);
%                 new_centers(i, :) = PopDec(randnum, :);
%                 rol_num = 0;
%             else
%                 new_centers(i, :) = mean(PopDec(center_index, :), 1);
%             end
%         end
%         if all(all(new_centers == centers))
%             break;
%         end
%         rol_num = rol_num+1;
%         centers = new_centers;
%     end
end