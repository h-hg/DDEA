function output = sort_dis(DataBase)

global Dim
[NP_DataBase, ~] = size(DataBase);

difference = DataBase(:, 1:Dim) - repmat(DataBase(1, 1:Dim), NP_DataBase, 1);
square = difference.^2;
summation = sum(square, 2);
distance = sqrt(summation);
DataBase = [DataBase, distance];
DataBase = sortrows(DataBase, Dim + 2);
DataBase = DataBase(:, 1:Dim + 1);
output = DataBase;