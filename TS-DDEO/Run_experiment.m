%------------------------------------------------------------------------
% This code is part of the program that produces the results in the following paper:
% Huixiang Zhen, Wenyin Gong, Ling Wang, Fei Ming, and Zuowen Liao. "Two-stage Data-driven Evolutionary Optimization for High-dimensional Expensive Problems", IEEE Transactions on Cybernetics, accepted, 2021.
% You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%----------------------------------------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     this code is used to run experiments on testset     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
addpath(genpath(pwd));

% Experiment Parameter
TestFuns = {  'ELLIPSOID'; 'ROSENBROCK'; 'ACKLEY'; 'GRIEWANK'; 'CEC05_f10'; 'CEC05_f16'; 'CEC05_f19'};    % The objective functions to be tested
dims = [30 50 100];                     % Dimensions to be tested
Runs = 1;                               % Number of runs

d = size(dims,2);
o = length(TestFuns);

bag_gsamp1 = {};                        % pack gsamp1
bag_time_cost = {};                     % pack time cost

% runs according to dims and objs.
for i = 1:d
    for j = 1:o
        fname = cell2mat(TestFuns(j));                  
        FUN=@(x) feval(fname,x); 
        [Xmin, Xmax] = variable_domain(fname); 
        LB = repmat((Xmin),1,dims(i));             
        UB = repmat((Xmax),1,dims(i),1);
        [ gsamp1,time_cost] = RUN_TSDDEO(Runs,dims(i),FUN, LB, UB, fname);
        bag_result(j,i) = {gsamp1(end)};
        bag_gsamp1(j,i) = {gsamp1};
        bag_time_cost(j,i) = {time_cost};
    end
end
save Result                     