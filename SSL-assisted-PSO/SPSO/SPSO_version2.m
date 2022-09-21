function [G_Optimal Iter_OptimalResult Time] = SPSO_original(func_id,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations)
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

if (func_id == 9 || func_id == 11 || func_id == 13)
    flag1 = 10;
else if (func_id == 10 || func_id == 12 || func_id == 14)
        flag1 = 19;
    end
end

for r=1:run_times
    initial_time=cputime;
    current_velocity = zeros(popsize,dimension);
    current_position=zeros(popsize,dimension);
    previous_position=zeros(popsize,dimension);
    current_fitness=zeros(popsize,1);
    pbest = zeros(popsize,dimension);
    pbest_fitness = zeros(popsize,1);
    gbest=zeros(1,dimension);
    gbest_fitness=zeros(1,1);
    previous_gbest_fitness=zeros(1,1);
    pbest_kept = zeros(popsize,dimension);
    pbestfitness_kept = zeros(popsize,1);
    
    gbest_output=[];
    evaltimes_output=[];
    archive_position = zeros(max_evaluations,dimension);
    archive_fitness = zeros(max_evaluations,1);
    n_archive = 1;

    evaltimes=0;
    iter = 1;
    
    n_nochange = 0;
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
            current_fitness(i,1) = benchmark_func(temp_position,flag1);
        else
            current_fitness(i,1)=fitness(temp_position,func_id);
        end
        pbest(i,:) = current_position(i,:);
        pbest_fitness(i,1) = current_fitness(i,1);
        archive_position(n_archive,:) = current_position(i,:);
        archive_fitness(n_archive,1) = current_fitness(i,1);
        n_archive = n_archive + 1;
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
     previous_gbest_fitness(1,1) = gbest_fitness(1,1);
     
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
%              for d=1:dimension
%                  avg(1,d) = 0;
%                  for k=1:n_archive-1
%                      avg(1,d) = avg(1,d) + archive_position(k,d);
%                  end
%                  avg(1,d) = avg(1,d)/popsize;
%              end
            
             for d=1:dimension
                 %current_position(i,d) = rand()*current_position(i,d) +  rand()*(pbest(selected_examplar,d) - current_position(i,d)) + rand()*(avg(1,d) - current_position(i,d));
                 previous_position(i,d) = current_position(i,d);
                 current_position(i,d) = rand()*current_position(i,d) +  rand()*(pbest(selected_examplar,d) - current_position(i,d));
%                  current_velocity(i,d) = 0.7293*(current_velocity(i,d) + 2.05*rand()*(pbest(i,d) - current_position(i,d)) + 2.05*rand()*(gbest(1,d) - current_position(i,d)));
%                  if current_velocity(i,d) < vmin
%                      current_velocity(i,d) = vmin;
%                  end
%                  if current_velocity(i,d) > vmax
%                      current_velocity(i,d) = vmax;
%                  end
%                  current_position(i,d) = current_position(i,d) + current_velocity(i,d);
                 if current_position(i,d) < xmin
                     current_position(i,d) = xmin;
                 end
                 if current_position(i,d) > xmax
                     current_position(i,d) = xmax;
                 end
             end
         end
         
         if n_nochange == 5
             for i=1:popsize
                 %real evaluation
                 temp_position=current_position(i,:);
                 %current_fitness(i,1)=fitness(temp_position,4);
                 if func_id >= 9 && func_id <= 14
                     current_fitness(i,1) = benchmark_func(temp_position,flag1);
                 else
                     current_fitness(i,1)=fitness(temp_position,func_id);
                 end
                 evaltimes=evaltimes+1;
                 gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
                 if current_fitness(i,1) < pbest_fitness(i,1)
                     pbest(i,:) = current_position(i,:);
                     pbest_fitness(i,1) = current_fitness(i,1);
                 end
                 for k=1:n_archive-1
                     dist_to_archive(k,1) = 0;
                     for d=1:dimension
                         dist_to_archive(k,1) = dist_to_archive(k,1) + (current_position(i,d) - archive_position(k,d))^2;
                     end
                     dist_to_archive(k,1) = sqrt(dist_to_archive(k,1));
                 end
                 [sort_value sort_id] = sort(dist_to_archive,'ascend');
                 if (sort_value(1,1) ~= 0)
                     archive_position(n_archive,:) = current_position(i,:);
                     archive_fitness(n_archive,1) = current_fitness(i,1);
                     n_archive = n_archive + 1;
                 end
             end
             n_nochange = 0;
         else
             for i=1:popsize
                 neighbor_threshold(i,1) = 0;
                 for d=1:dimension
                     neighbor_threshold(i,1) = neighbor_threshold(i,1) + (current_position(i,d) - previous_position(i,d))^2;
                 end
                 neighbor_threshold(i,1) = sqrt(neighbor_threshold(i,1));
             end
             
             n_currarchive = n_archive-1;
             dist_archive = [];
             for i=1:popsize
                 for j=1:n_currarchive
                     dist_archive(i,j) = 0;
                     for d=1:dimension
                         dist_archive(i,j) = dist_archive(i,j) + (current_position(i,d) - archive_position(j,d))^2;
                     end
                     dist_archive(i,j) = dist_archive(i,j)^0.5;
                 end
             end
             for i=1:popsize
                 [sort_value sort_id] = sort(dist_archive(i,j),'ascend');
                 min_dist = sort_value(1,1);
                 if min_dist < 1e-3
                     if archive_fitness(sort_id(1,1),1) < pbest_fitness(i,1)
                         %real evaluation
                         temp_position=current_position(i,:);
                         %current_fitness(i,1)=fitness(temp_position,4);
                         if func_id >= 9 && func_id <= 14
                             current_fitness(i,1) = benchmark_func(temp_position,flag1);
                         else
                             current_fitness(i,1)=fitness(temp_position,func_id);
                         end
                         evaltimes=evaltimes+1;
                         gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
                         if current_fitness(i,1) < pbest_fitness(i,1)
                             pbest(i,:) = current_position(i,:);
                             pbest_fitness(i,1) = current_fitness(i,1);
                         end
                         if abs(current_fitness(i,1) - archive_fitness(sort_id(1,1),1)) > 1e-3
                             archive_position(n_archive,:) = current_position(i,:);
                             archive_fitness(n_archive,1) = current_fitness(i,1);
                             n_archive = n_archive + 1;
                         end
                     end
                 else
                     [close_dist close_id] = find(dist_archive(i,:) <= neighbor_threshold(i,1));
                     n_close_id = size(close_id',1);
                     [better_value better_id] = find(archive_fitness(close_id(1,:),1) <= pbest_fitness(i,1));
                     if size(better_id,1) >= fix(n_close_id/3)+1
                         %real evaluation
                         temp_position=current_position(i,:);
                         %current_fitness(i,1)=fitness(temp_position,4);
                         if func_id >= 9 && func_id <= 14
                             current_fitness(i,1) = benchmark_func(temp_position,flag1);
                         else
                             current_fitness(i,1)=fitness(temp_position,func_id);
                         end
                         evaltimes=evaltimes+1;
                         gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
                         if current_fitness(i,1) < pbest_fitness(i,1)
                             pbest(i,:) = current_position(i,:);
                             pbest_fitness(i,1) = current_fitness(i,1);
                         end
                         archive_position(n_archive,:) = current_position(i,:);
                         archive_fitness(n_archive,1) = current_fitness(i,1);
                         n_archive = n_archive + 1;
                     end
                 end
             end
         end
         for i=1:popsize
             if pbest_fitness(i,1)<gbest_fitness(1,1)
                 gbest(1,:) = pbest(i,:);
                 gbest_fitness(1,1) = pbest_fitness(i,1);
             end
         end
         
         if gbest_fitness(1,1) == previous_gbest_fitness(1,1)
             n_nochange = n_nochange + 1;
         else
             n_nochange = 0;
             previous_gbest_fitness(1,1) = gbest_fitness(1,1);
         end
         
         iter = iter + 1
         gbest_fitness
         evaltimes
         n_archive
     end
     Time(r,1) = cputime - initial_time;
     G_Optimal(r,:) = gbest_fitness(1,1);
     for k=1:max_evaluations-popsize
         Iter_OptimalResult(r,k) = gbest_output(1,k);
     end
end
         
end