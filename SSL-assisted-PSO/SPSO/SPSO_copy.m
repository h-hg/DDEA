function [G_Optimal Iter_OptimalResult Time] = SPSO_w(func_id,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations)
% clc
% clear all
% global dimension
% global popsize
global max_length
global evaltimes

global initial_flag

initial_flag=0;
n_fiteval=0; 
format short; 
max_iteration=zeros(run_times,1);

N_FitEval = [];
Iter_OptimalResult = [];
G_Optimal = [];
Time = [];

if (func_id == 9 || func_id == 11 || func_id == 13)
    flag1 = 10;
else if (func_id == 10 || func_id == 12 || func_id == 14)
        flag1 = 9;
    end
end

c1 = 2.05;
c2 = 2.05;
w = 0.7;
max_node=dimension;
max_evaluations = 11*dimension;

for r=1:run_times
    initial_time=cputime;
    current_velocity = zeros(popsize,dimension);
    current_position=zeros(popsize,dimension);
    current_fitness=zeros(popsize,1);
    pbest = zeros(popsize,dimension);
    pbest_fitness = zeros(popsize,1);
    gbest=zeros(1,dimension);
    gbest_fitness=zeros(1,1);
    
    evaluated = ones(popsize,1);
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
    
     while evaltimes<max_evaluations
         [sort_pbestfit sort_pbestid] = sort(pbest_fitness);
         %[sort_pbestfit sort_pbestid] = sort(current_fitness);
         for i=1:popsize
             sortid = find(sort_pbestid == i);
             for d=1:dimension
                 randset = randperm(sortid);
                 selected_pbest = sort_pbestid(randset(1,1));
                 current_velocity(i,d) = rand*current_velocity(i,d) + rand*(pbest(selected_pbest,d)-current_position(i,d));
                 %current_velocity(i,d) = rand*current_velocity(i,d) + rand*(current_position(selected_pbest,d)-current_position(i,d));
                 %current_velocity(i,d) = 0.729*(current_velocity(i,d) + c1*rand*(pbest(i,d)-current_position(i,d)) + c2*rand*(gbest(1,d)-current_position(i,d)));
                 if current_velocity(i,d) < vmin
                     current_velocity(i,d) = vmin;
                 end
                 if current_velocity(i,d) > vmax
                     current_velocity(i,d) = vmax;
                 end
                 current_position(i,d) = current_position(i,d) + current_velocity(i,d);
                 if current_position(i,d) < xmin
                     current_position(i,d) = xmin;
                     %current_position(i,d) = rand*(xmax-xmin)+xmin;
                 end
                 if current_position(i,d) > xmax
                     current_position(i,d) = xmax;
                     %current_position(i,d) = rand*(xmax-xmin)+xmin;
                 end
             end
         end
         
         for d=1:dimension
             upper_d(1,d) = max(current_position(:,d));
             %upper_dim(1,d) = min(archive_position(find(archive_position(:,d)>upper_d(1,d)),d));
             lower_d(1,d) = min(current_position(:,d));
             %lower_dim(1,d) = max(archive_position(find(archive_position(:,d)<lower_d(1,d)),d));
         end
         trainset = [];
         trainset_fitness = [];
         train_num = 1;
         for arch=1:n_archive-1
             flag = 1;
             for d=1:dimension
                 if archive_position(arch,d)<lower_d(1,d) || archive_position(arch,d)>upper_d(1,d)
                     flag = 0;
                     break;
                 end
             end
             if flag == 1
                 if size(find(abs(trainset_fitness - archive_fitness(arch,1)) < 1e-3),1) == 0
                     trainset(train_num,:) = archive_position(arch,:);
                     trainset_fitness(train_num,1) = archive_fitness(arch,1);
                     train_num = train_num+1;
                 else
                     same = 1;
                     for t=1:train_num-1
                         if size(find(abs(trainset(t,:)-archive_position(arch,:))>1e-1),1)>1
                             same = 0;
                             break;
                         end
                     end
                     if same == 0
                         trainset(train_num,:) = archive_position(arch,:);
                         trainset_fitness(train_num,1) = archive_fitness(arch,1);
                         train_num = train_num+1;
                     end
                 end
             end
         end
         train_num = train_num-1;
         
         if train_num<dimension+1
             trainset = [];
             trainset_fitness = [];
             train_num = 1;
             for arch=1:n_archive-1
                 if size(find(abs(trainset_fitness - archive_fitness(arch,1)) < 1e-3),1) == 0
                     trainset(train_num,:) = archive_position(arch,:);
                     trainset_fitness(train_num,1) = archive_fitness(arch,1);
                     train_num = train_num+1;
                 else
                     same = 1;
                     for t=1:train_num-1
                         if size(find(abs(trainset(t,:)-archive_position(arch,:))>1e-1),1)>1
                             same = 0;
                             break;
                         end
                     end
                     if same == 0
                         trainset(train_num,:) = archive_position(arch,:);
                         trainset_fitness(train_num,1) = archive_fitness(arch,1);
                         train_num = train_num+1;
                     end
                 end
             end
             train_num = train_num-1;
         end
         
         for d=1:dimension
%              max_range(d,1)=max(trainset(:,d));
%              min_range(d,1)=min(trainset(:,d));
             max_range(d,1)=max(current_position(:,d));
             min_range(d,1)=min(current_position(:,d));
         end
         %     spread=max(max_range-min_range);
         spread=0;
         for d=1:dimension
             spread=spread+(max_range(d,1)-min_range(d,1))^2;
         end
         spread=sqrt(spread);
         %trainmodel=newrb(trainset',trainset_fitness',0.01,spread,max_node,1);
         trainmodel = newrbe(trainset',trainset_fitness',spread);
         
         n_currarchive = n_archive-1;
         min_diff = abs(max(archive_fitness)-min(archive_fitness));
         max_diff = 0;
         mindiff_id = 0;
         maxdiff_id = 0;
         min_neighbordist = 0;
         min_neighid = 0;         
         for i=1:popsize
             current_fitness(i,:) = sim(trainmodel,current_position(i,:)');
             if current_fitness(i,:) < pbest_fitness(i,1)
                 %current_fitness(i,:) = approx_fitness;
                 pbest(i,:) = current_position(i,:);
                 pbest_fitness(i,1) = current_fitness(i,1);
                 evaluated(i,1) = 0;
                 dist=zeros(n_archive-1,1);
                 for k=1:n_archive-1
                     temp_dist = 0;
                     for d=1:dimension
                         temp_dist = temp_dist + (pbest(i,d)-archive_position(k,d))^2;
                     end
                     dist(k,1) = sqrt(temp_dist);
                 end
                 [min_distvalue min_idvalue] = min(dist);
                 if min_distvalue>min_neighbordist
                     min_neighbordist = min_distvalue;
                     min_neighborid = min_neighid;
                 end
%                  archive_position(n_archive,:) = current_position(i,:);
%                  archive_fitness(n_archive,1) = current_fitness(i,1);
%                  n_archive = n_archive + 1;
             else if current_fitness(i,:) > pbest_fitness(i,1)
                     %current_fitness(i,1) = approx_fitness;
                     if abs(gbest_fitness-current_fitness(i,1))<min_diff
                         min_diff = abs(gbest_fitness-current_fitness(i,1));
                         mindiff_id = i;
                     end
                     if abs(gbest_fitness-current_fitness(i,1))>max_diff
                         max_diff = abs(gbest_fitness-current_fitness(i,1));
                         maxdiff_id = i;
                     end
                 end
             end
         end
%          if maxdiff_id ~= 0
%              if func_id >= 9 && func_id <= 14
%                  %current_fitness(mindiff_id,1) = benchmark_func(current_position(mindiff_id,:),flag1);
%                  current_fitness(maxdiff_id,1) = benchmark_func(current_position(maxdiff_id,:),flag1);
%              else
%                  %current_fitness(mindiff_id,1)=fitness(current_position(mindiff_id,:),func_id);
%                  current_fitness(maxdiff_id,1) = fitness(current_position(maxdiff_id,:),func_id);
%              end
%              evaltimes=evaltimes+1;
%              gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
% %              archive_position(n_archive,:) = current_position(mindiff_id,:);
% %              archive_fitness(n_archive,1) = current_fitness(mindiff_id,1);
% %              n_archive = n_archive + 1;
%              archive_position(n_archive,:) = current_position(maxdiff_id,:);
%              archive_fitness(n_archive,1) = current_fitness(maxdiff_id,1);
%              n_archive = n_archive + 1;
% %              if current_fitness(mindiff_id,1) < pbest_fitness(mindiff_id,1)
% %                  pbest(mindiff_id,:) = current_position(mindiff_id,:);
% %                  pbest_fitness(mindiff_id,1) = current_fitness(mindiff_id,1);
% %                  evaluated(mindiff_id,1) = 1;
% %              end
%              if current_fitness(maxdiff_id,1) < pbest_fitness(maxdiff_id,1)
%                  pbest(maxdiff_id,:) = current_position(maxdiff_id,:);
%                  pbest_fitness(maxdiff_id,1) = current_fitness(maxdiff_id,1);
%                  evaluated(maxdiff_id,1) = 1;
%              end
%              
% %              if min_neighid ~= 0 & maxdiff_id ~= min_neighid
% %                  archive_position(n_archive,:) = current_position(min_neighid,:);
% %                  archive_fitness(n_archive,1) = current_fitness(min_neighid,1);
% %                  n_archive = n_archive + 1;
% %              end
%          end
         if mindiff_id ~= 0
             if func_id >= 9 && func_id <= 14
                 current_fitness(mindiff_id,1) = benchmark_func(current_position(mindiff_id,:),flag1);
                 %current_fitness(maxdiff_id,1) = benchmark_func(current_position(maxdiff_id,:),flag1);
             else
                 current_fitness(mindiff_id,1)=fitness(current_position(mindiff_id,:),func_id);
                 %current_fitness(maxdiff_id,1) = fitness(current_position(maxdiff_id,:),func_id);
             end
             evaltimes=evaltimes+1;
             gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
%              archive_position(n_archive,:) = current_position(mindiff_id,:);
%              archive_fitness(n_archive,1) = current_fitness(mindiff_id,1);
%              n_archive = n_archive + 1;
             archive_position(n_archive,:) = current_position(mindiff_id,:);
             archive_fitness(n_archive,1) = current_fitness(mindiff_id,1);
             n_archive = n_archive + 1;
%              if current_fitness(mindiff_id,1) < pbest_fitness(mindiff_id,1)
%                  pbest(mindiff_id,:) = current_position(mindiff_id,:);
%                  pbest_fitness(mindiff_id,1) = current_fitness(mindiff_id,1);
%                  evaluated(mindiff_id,1) = 1;
%              end
             if current_fitness(mindiff_id,1) < pbest_fitness(mindiff_id,1)
                 pbest(mindiff_id,:) = current_position(mindiff_id,:);
                 pbest_fitness(mindiff_id,1) = current_fitness(mindiff_id,1);
                 evaluated(mindiff_id,1) = 1;
             end
         end
         
         if min_neighborid ~= 0 %& maxdiff_id ~= min_neighid
             %if min_neighid ~= 0 & mindiff_id ~= min_neighid & maxdiff_id ~= min_neighid
             archive_position(n_archive,:) = current_position(min_neighborid,:);
             archive_fitness(n_archive,1) = current_fitness(min_neighborid,1);
             n_archive = n_archive + 1;
         end
         
         [minvalue minid] = min(pbest_fitness);
         if minid == min_neighid
             n_archive = n_archive-1;   
         end
         if minvalue<gbest_fitness(1,1) 
             if evaluated(minid,1)==0
                 if func_id >= 9 && func_id <= 14
                     pbest_fitness(minid,1) = benchmark_func(pbest(minid,:),flag1);
                 else
                     pbest_fitness(minid,1)=fitness(pbest(minid,:),func_id);
                 end
                 evaluated(minid,1) = 1;
                 current_fitness(minid,1) = pbest_fitness(minid,1);
                 archive_position(n_archive,:) = pbest(minid,:);
                 archive_fitness(n_archive,1) = pbest_fitness(minid,:);
                 n_archive = n_archive + 1;
                 evaltimes=evaltimes+1;
                 if pbest_fitness(minid,1) < gbest_fitness
                     gbest(1,:) = pbest(minid,:);
                     gbest_fitness(1,1) = pbest_fitness(minid,1);
                 end
             else
                 gbest(1,:) = pbest(minid,:);
                 gbest_fitness(1,1) = pbest_fitness(minid,1);
             end
         end
%          for i=1:popsize
%              if pbest_fitness(i,1)<gbest_fitness(1,1)
%                  gbest(1,:) = pbest(i,:);
%                  gbest_fitness(1,1) = pbest_fitness(i,1);
%              end
%          end
%          trainset = archive_position;
%          trainset_fitness = archive_fitness;
%          for d=1:dimension
%              temp_trainspread(d,1) = max(current_position(:,d))-min(current_position(:,d));
%          end
%          trainspread = mean(temp_trainspread)*0.1;
%          trainmodel = newrbe(trainset',trainset_fitness',trainspread);
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