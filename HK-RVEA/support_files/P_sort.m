function [FrontValue,MaxFront] = P_sort(FunctionValue,operation)
% 进行非支配排序
% 输入: FunctionValue, 待排序的种群(目标空间)
%       Operation,     可指定仅排序第一个面,排序前一半个体,或是排序所有的个体, 默认为排序所有的个体
% 输出: FrontValue, 排序后的每个个体所在的前沿面编号, 未排序的个体前沿面编号为inf
%       MaxFront,   排序的最大前面编号

    if nargin<2
        kinds=1;
    elseif strcmp(operation,'half')
        kinds=2;
    elseif strcmp(operation,'first')
        kinds=3;
    else
        kinds=1;
    end
    
    [N,M] = size(FunctionValue);
    MaxFront = 0;
    Sorted = false(1,N);
    [FunctionValue,rank] = sortrows(FunctionValue);
    FrontValue = zeros(1,N) + inf;
    while (kinds==1 && sum(Sorted)<N) || (kinds==2 && sum(Sorted)<N/2) || (kinds==3 && MaxFront<1)
        MaxFront = MaxFront + 1;
        ThisFront = false(1,N);
        for i = 1 : N
            if ~Sorted(i)
                x = 0;
                for j = 1 : N
                    if ThisFront(j)
                        x = 2;
                        for j2 = 2 : M
                            if FunctionValue(i,j2) < FunctionValue(j,j2)
                                x = 0;
                                break;
                            end
                        end
                        if x == 2
                            break;
                        end
                    end
                end
                if x ~= 2
                    ThisFront(i) = true;
                    Sorted(i) = true;
                end
            end
        end
        FrontValue(rank(ThisFront)) = MaxFront;
    end
end