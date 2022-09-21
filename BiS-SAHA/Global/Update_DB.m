function Update_DB(eval_individual)

global Dim
global DataBase original_DB
global num_eval
global best_solution
global grid
if ~isempty(setdiff(eval_individual(:, 1:Dim), DataBase(:, 1:Dim), 'rows'))
    num_eval = num_eval + 1;
    eval_individual = finess(eval_individual(:, 1:Dim));
    DataBase = [DataBase; eval_individual];
    original_DB = [eval_individual; original_DB];
    DataBase = sortrows(DataBase, Dim + 1);
    best_solution = [best_solution; DataBase(1, :)];
end