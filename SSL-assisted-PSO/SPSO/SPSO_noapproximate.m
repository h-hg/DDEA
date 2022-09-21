function [G_Optimal Iter_OptimalResult Time] = SPSO_noapproximate(func_id,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations)
% clc
% clear all
% global dimension
% global popsize
global max_length
global evaltimes

global initial_flag

initial_flag=0;
% dimension=50;
% popsize=30;
%run_times=20;
% max_evaluations=1000;
n_fiteval=0; 
%arch_size = 5*dimension;
format short; 
max_iteration=zeros(run_times,1);

N_FitEval = [];
Iter_OptimalResult = [];
G_Optimal = [];
Time = [];

if func_id == 9 || func_id == 11 || func_id == 13
    flag = 10;
else if func_id == 10 || func_id == 12 || func_id == 14
        flag = 9;
    end
end

for r=1:run_times
    initial_time=cputime;
    current_velocity = zeros(popsize,dimension);
    current_position=zeros(popsize,dimension);
    current_fitness=zeros(popsize,1);
    pbest = zeros(popsize,dimension);
    pbest_fitness = zeros(popsize,1);
    gbest=zeros(1,dimension);
    gbest_fitness=zeros(1,1);
    
    pbest_kept = zeros(popsize,dimension);
    pbestfitness_kept = zeros(popsize,1);
    
    gbest_output=[];
    evaltimes_output=[];
    archive_position = zeros(max_evaluations,dimension);
    archive_fitness = zeros(max_evaluations,1);
    n_archive = 1;

    evaltimes=0;
    iter = 1;
    for i=1:popsize
        for d=1:dimension
            current_velocity(i,d) = rand*(vmax - vmin)+vmin;
        end
    end
    %generate initial population using Latin Hypervolume method
    [current_position] = latinhv(popsize,dimension,xmin,xmax);

    for i=1:popsize
        temp_position=current_position(i,:);
        if func_id >= 9 && func_id <= 14
            current_fitness(i,1) = benchmark_func(temp_position,flag);
        else
            current_fitness(i,1)=fitness(temp_position,func_id);
        end
        pbest(i,:) = current_position(i,:);
        pbest_fitness(i,1) = current_fitness(i,1);
        evaltimes=evaltimes+1;
    end

    evaltimes_output(1)=evaltimes;

     gbest(1,:)=current_position(1,:);
     gbest_fitness(1,1)=current_fitness(1,1);
     for i=2:popsize
         if current_fitness(i,1)<gbest_fitness(1,1)
             gbest(1,:)=current_position(i,:);
             gbest_fitness(1,1)=current_fitness(i,1);
         end
     end
     gbest_output(1)=gbest_fitness(1,1); 
     
     while evaltimes<1000
         [sort_value sort_id] = sort(pbest_fitness);
         for i=1:popsize
             %better_id = find(sort_value<mean(pbest_fitness));
             %better_id = find(sort_value<(max(pbest_fitness)+min(pbest_fitness))/2);
             %selected_id = fix(rand()*size(better_id,1));
             selected_id = fix(rand()*(popsize/3));
             while selected_id == 0
                 %selected_id = fix(rand()*size(better_id,1));
                 selected_id = fix(rand()*(popsize/3));
             end
             selected_examplar = sort_id(selected_id);
             same_position = 1;
             while same_position == 1
                for d=1:dimension
                    if current_position(i,d) ~= pbest(selected_examplar,d)
                        same_position = 0;
                        break;
                    end
                end
                if same_position == 1
                    %better_id = find(sort_value<mean(pbest_fitness));
                    %better_id = find(sort_value<(max(pbest_fitness)+min(pbest_fitness))/2);
                    %selected_id = fix(rand()*size(better_id,1));
                    selected_id = fix(rand()*(popsize/3));
                    while selected_id == 0
                        %selected_id = fix(rand()*size(better_id,1));
                        selected_id = fix(rand()*(popsize/3));
                    end
                    selected_examplar = sort_id(selected_id);
                end
             end
             for d=1:dimension
                 %current_position(i,d) = rand()*current_position(i,d) +  rand()*(pbest(selected_examplar,d) - current_position(i,d)) + rand()*(avg(1,d) - current_position(i,d));
                 current_position(i,d) = rand()*current_position(i,d) +  rand()*(pbest(selected_examplar,d) - current_position(i,d));
                 if current_position(i,d) < xmin
                     current_position(i,d) = xmin;
                 end
                 if current_position(i,d) > xmax
                     current_position(i,d) = xmax;
                 end
             end
         end
         
         for i=1:popsize
             temp_position=current_position(i,:);
             %current_fitness(i,1)=fitness(temp_position,4);
             if func_id >= 9 && func_id <= 14
                 current_fitness(i,1) = benchmark_func(temp_position,flag);
             else
                 current_fitness(i,1)=fitness(temp_position,func_id);
             end
             evaltimes=evaltimes+1;
             gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
             if current_fitness(i,1) < pbest_fitness(i,1)
                 pbest(i,:) = current_position(i,:);
                 pbest_fitness(i,1) = current_fitness(i,1);
             end
         end

         for i=1:popsize
             if pbest_fitness(i,1)<gbest_fitness(1,1)
                 gbest(1,:) = pbest(i,:);
                 gbest_fitness(1,1) = pbest_fitness(i,1);
             end
         end
         iter = iter + 1
         gbest_fitness
         evaltimes
     end
     Time(r,1) = cputime - initial_time;
     G_Optimal(r,:) = gbest_fitness(1,1);
     for k=1:max_evaluations-popsize
         Iter_OptimalResult(r,k) = gbest_output(1,k);
     end
end
         
end