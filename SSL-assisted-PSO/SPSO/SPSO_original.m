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
         
%          min_position = zeros(1,dimension);
%          max_position = zeros(1,dimension);
%          for d=1:dimension
%              min_position(1,d) = min(current_position(:,d));
%              max_position(1,d) = max(current_position(:,d));
%          end
%          max_dist = 0;
%          for d=1:dimension
%              max_dist = max_dist + (max_position(1,d) - min_position(1,d))^2;
%          end
%          max_dist = max_dist^0.5;
%          neighbor_threshold = max_dist/2;
%          dist_archive = zeros(popsize,n_archive - 1);
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
         closest_value = [];
         closest_ind = [];
         for j=1:n_currarchive
             [closest_value(j,1) closest_ind(j,1)] = min(dist_archive(:,j));
         end
         for i=1:popsize
             n_closest = [];
             n_closest = find(closest_ind == i);
             if size(n_closest,1)>=3
             [min_value minid] = min(dist_archive(i,n_closest));
             neighbor_id = find(dist_archive(i,:) <= min_value);
             better_inds = find(archive_fitness(neighbor_id,1)<pbest_fitness(i,1));
             if size(better_inds,1) > (1+fix(size(neighbor_id,1))/2); 
                 %real evaluation
                 temp_position=current_position(i,:);
                     %current_fitness(i,1)=fitness(temp_position,4);
                     if func_id >= 9 && func_id <= 14
                         current_fitness(i,1) = benchmark_func(temp_position,flag1);
                     else
                         current_fitness(i,1)=fitness(temp_position,func_id);
                     end
                     same_index = find(abs(current_fitness(i,1) - archive_fitness(1:n_archive-1,1))<1e-3);
                     for k=1:size(same_index,1)
                         archive_indid = same_index(k,1);
                         dist_archive1 = 0;
                         for d=1:dimension
                             dist_archive1 = dist_archive1 + (current_position(i,d) - archive_position(archive_indid,d))^2;
                         end
                         dist_archive1 = sqrt(dist_archive1);
                         if dist_archive1 > 1e-3
                             archive_position(n_archive,:) = current_position(i,:);
                             archive_fitness(n_archive,1) = current_fitness(i,1);
                             n_archive = n_archive + 1;
                         else
                             if current_fitness(i,1) < archive_fitness(archive_indid,1)
                                 archive_fitness(archive_indid,1) = current_fitness(i,1);
                                 archive_position(archive_indid,:) = current_position(i,:);
                             end
                         end
                     end
                     evaltimes=evaltimes+1;
                     gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
                     if current_fitness(i,1) < pbest_fitness(i,1)
                         pbest(i,:) = current_position(i,:);
                         pbest_fitness(i,1) = current_fitness(i,1);
                     end
             end
             else
                 temp_position=current_position(i,:);
                     %current_fitness(i,1)=fitness(temp_position,4);
                     if func_id >= 9 && func_id <= 14
                         current_fitness(i,1) = benchmark_func(temp_position,flag1);
                     else
                         current_fitness(i,1)=fitness(temp_position,func_id);
                     end
                     same_index = find(abs(current_fitness(i,1) - archive_fitness(1:n_archive-1,1))<1e-3);
                     for k=1:size(same_index,1)
                         archive_indid = same_index(k,1);
                         dist_archive1 = 0;
                         for d=1:dimension
                             dist_archive1 = dist_archive1 + (current_position(i,d) - archive_position(archive_indid,d))^2;
                         end
                         dist_archive1 = sqrt(dist_archive1);
                         if dist_archive1 > 1e-3
                             archive_position(n_archive,:) = current_position(i,:);
                             archive_fitness(n_archive,1) = current_fitness(i,1);
                             n_archive = n_archive + 1;
                         else
                             if current_fitness(i,1) < archive_fitness(archive_indid,1)
                                 archive_fitness(archive_indid,1) = current_fitness(i,1);
                                 archive_position(archive_indid,:) = current_position(i,:);
                             end
                         end
                     end
                     evaltimes=evaltimes+1;
                     gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
                     if current_fitness(i,1) < pbest_fitness(i,1)
                         pbest(i,:) = current_position(i,:);
                         pbest_fitness(i,1) = current_fitness(i,1);
                     end
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