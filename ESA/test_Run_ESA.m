% Experiments for Paper
clc;
clear;
warning('off');
addpath(genpath(pwd));
TestFuns = {'CEC05_F10';'CEC05_F19'; 'ACKLEY';  'GRIEWANK'; 'ROSENBROCK'; 'ELLIPSOID'; 'CEC05_F16';};    % The objective functions to be tested
dims = [30  100 50 ];                       % Dimensions to be tested
Runs = 1;                                  % Number of runs

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
        [ gsamp1,time_cost] = RUN_ESA(Runs,dims(i),FUN, LB, UB, fname);
        
        % Each line contains result of a function with different dimensions
        bag_result(j,i) = {gsamp1(end)};
        bag_gsamp1(j,i) = {gsamp1};
        bag_time_cost(j,i) = {time_cost};
    end
end
save Result                     