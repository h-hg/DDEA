% A Kriging-assisted evolutionary algorithm based on incremental learning
%
% Dawei Zhan and Huanlai Xing. A Fast Kriging-Assisted Evolutionary Algorithm
% Based on Incremental Learning.IEEE Transactions on Evolutionary Computation,
% doi: 10.1109/TEVC.2021.3067015.
% the implementation of the incremental Kriging was based on the code of [1]
% [1] A. I. J. Forrester, A. Sobester and A. J. Keane (2008). Engineering design
% via surrogate modelling: a practical guide, John Wiley & Sons.
% 
% coded by Dawei Zhan (zhandawei at swjtu dot edu dot cn)
clearvars;clc;close all;
% settings of the test problem
fun_name = 'Ellipsoid';
num_vari = 100;
lower_bound = -5.12*ones(1,num_vari); 
upper_bound = 5.12*ones(1,num_vari);
% algorithm settings
num_initial = 100;
max_evaluation = 1000;
% DE settings
NP = 50;
CR = 0.8;
F = 0.8;
% number of current generation
generation = 1;
% generate random samples
sample_x = lhsdesign(num_initial, num_vari,'criterion','maximin','iterations',1000).*(upper_bound - lower_bound) + lower_bound;
sample_y = feval(fun_name, sample_x);
evaluation =  size(sample_x,1);
% best objectives in each generation
fmin_record = zeros(max_evaluation - evaluation + 1,1);
% the first DE population
[~,index] = sort(sample_y);
pop_vari = sample_x(index(1:NP),:);
pop_obj = sample_y(index(1:NP),:);
fmin = sample_y(index(1));
xmin = sample_x(index(1),:);
fmin_record(generation,:) = fmin;
time = 0;
% print the iteration information
fprintf('IKAEA on %d-D %s , generation: %d, evaluation: %d, best: %0.4g, total training time: %0.4g\n',num_vari,fun_name,generation,evaluation,fmin,time);
% the evoluation of the population
while evaluation < max_evaluation
    % build the Kriging model
    tic;
    if generation == 1
        % initial learning
        kriging_model = kriging_theta_train(sample_x,sample_y,lower_bound,upper_bound,1*ones(1,num_vari),0.000001*ones(1,num_vari),100*ones(1,num_vari));
    else
        % incremental learning
        kriging_model = kriging_incremental(kriging_model,infill_x,infill_y,sample_x(1:end-1,:),sample_y(1:end-1,:),lower_bound,upper_bound);
    end
    % mutation
    pop_mutation = zeros(NP,num_vari);
    for ii = 1 : NP
        avail_num = (1:NP);
        avail_num(ii) = [];
        % randomly generate three different integers
        r = avail_num(randperm(length(avail_num),2));
        pop_mutation(ii,:) = xmin + F*(pop_vari(r(1),:)-pop_vari(r(2),:));
        % check the bound constraints, randomly re-initialization
        if any(pop_mutation(ii,:)<lower_bound) || any(pop_mutation(ii,:)>upper_bound)
            pop_mutation(ii,:) = lower_bound + rand(1,num_vari).*(upper_bound-lower_bound);
        end
    end
    % crossover
    rand_matrix = rand(NP,num_vari);
    temp = randi(num_vari,NP,1);
    for ii = 1 : NP
        rand_matrix(ii,temp(ii)) = 0;
    end
    mui = rand_matrix < CR;
    pop_trial = pop_mutation.*mui + pop_vari.*(1-mui);
    % individuals too close to sampled points are removed
    pop_candi = [];
    for ii = 1 : NP
        if min(sqrt(sum((pop_trial(ii,:) - [sample_x;pop_candi]).^2,2)))>1E-6
            pop_candi = [pop_candi;pop_trial(ii,:)];
        end
    end
    % select infill samples (expected improvement criterion)
    [u,s] = kriging_predictor(pop_candi,kriging_model,sample_x,sample_y,lower_bound,upper_bound);
    [max_EI,ind] = max((fmin-u).*normcdf((fmin-u)./s)+s.*normpdf((fmin-u)./s));
    infill_x = pop_candi(ind,:);
    infill_y = feval(fun_name, infill_x);
    % replacement
    pop_vari_all = [pop_vari;infill_x];
    pop_obj_all = [pop_obj;infill_y];
    [~,index] = sort(pop_obj_all,'ascend');
    pop_vari = pop_vari_all(index(1:NP),:);
    pop_obj = pop_obj_all(index(1:NP),:);
    % update database
    sample_x = [sample_x;infill_x];
    sample_y = [sample_y;infill_y];
    % update the evaluation number and generation number
    generation = generation + 1;
    evaluation = evaluation + size(infill_x,1);
    % update current minimum
    [fmin,index] = min(sample_y);
    fmin_record(generation,:) = fmin;
    xmin = sample_x(index,:);
    time = time + toc;
    % print the iteration information
    fprintf('IKAEA on %d-D %s , generation: %d, evaluation: %d, best: %0.4g, total training time: %0.4g\n',num_vari,fun_name,generation,evaluation,fmin,time)
end


