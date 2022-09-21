%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     this code is used to run experiments on testset     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
addpath(genpath(pwd));

% Experiment Parameter
Runs = 1;                                                                   % Run number
dims = [30, 50, 100, 200];                                                  % different dimensions 
funarr = { 'Ellipsoid';};                                                   % ROSENBROCK ACKLEY GRIEWANK Ellipsoid CEC05_f10 CEC05_f16 CEC05_f19
o = size(funarr,1);                                                         % Number of objective function
d = size(dims,2);                                                           % Number of different dimensions 
bag_gsamp1 = {};
bag_time_cost = {};

% Runs according to dims and funarr.
for i = 1:d
    for j = 1:o
        fname = funarr{j};
        FUN=@(x) feval(fname,x);
        dim = dims(i);
        [Xmin, Xmax] = variable_domain(fname);
        LB = repmat((Xmin),1,dim);
        UB = repmat((Xmax),1,dim,1);
        
        % Run DESO
        [ gsamp1,time_cost] = RUN_DESO(FUN, Runs, dim, LB, UB, fname);
        
        % Each line contains the result of an objective function with different dimensions
        bag_result(j,i) = {gsamp1(end)};
        bag_gsamp1(j,i) = {gsamp1};
        bag_time_cost(j,i) = {time_cost};
    end
end
save Result