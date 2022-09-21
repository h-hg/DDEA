function output = init_pop_by_Clustering(DataBase)
% 2019-11-22 
%by:zhihai ren
global Size_pop Dim
num_class = 10;

T = clusterdata(DataBase(:, 1:Dim),'maxclust',num_class);
for i = 1:num_class
    class{i} = find(T == i);
end
sub_DataBase = [];
[NP_sub_DataBase, ~] = size(sub_DataBase);
while NP_sub_DataBase < Size_pop
    for i = 1:num_class
        current_class = class{i};
        NP_class = length(current_class);
        if NP_class > 0
            randIndex = randperm(NP_class);
            sub_DataBase = [sub_DataBase; DataBase(current_class(randIndex(1)), :)];
            current_class(randIndex(1)) = [];
            class{i} = current_class;
        end
        if NP_sub_DataBase >= Size_pop
            break;
        end
    end
    [NP_sub_DataBase, ~] = size(sub_DataBase);
end

output = sub_DataBase(1:Size_pop, :);