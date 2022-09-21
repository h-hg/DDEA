function [G_Optimal Iter_OptimalResult Time] = SPSO_w(func_id,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations)
%function [G_Optimal Time] = SPSO_w(func_id,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations)
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
max_node=8;
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
    archive_position = [];%zeros(max_evaluations,dimension);
    archive_fitness = [];%zeros(max_evaluations,1);
    label_position = [];
    label_fitness = [];
    unlabel_position = [];
    unlabel_fitness = [];
    n_unlabel = 0;
    n_label = 0;
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
        label_position(n_label+1,:) = current_position(i,:);
        label_fitness(n_label+1,1) = current_fitness(i,:);
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
         avg = mean(current_position);
         for i=1:popsize
             sortid = find(sort_pbestid == i);
             for d=1:dimension
                 randset = randperm(sortid);
                 selected_pbest = sort_pbestid(randset(1,1));
                 current_velocity(i,d) = rand*current_velocity(i,d) + rand*(pbest(selected_pbest,d)-current_position(i,d))+0.001*rand*(avg(1,d)-current_position(i,d));
                 if current_velocity(i,d) < vmin
                     current_velocity(i,d) = vmin;
                 end
                 if current_velocity(i,d) > vmax
                     current_velocity(i,d) = vmax;
                 end
                 current_position(i,d) = current_position(i,d) + current_velocity(i,d);
                 if current_position(i,d) < xmin
                     current_position(i,d) = xmin;
                 end
                 if current_position(i,d) > xmax
                     current_position(i,d) = xmax;
                 end
             end
         end
         
         for d=1:dimension
             upper_d(1,d) = max(current_position(:,d));
             lower_d(1,d) = min(current_position(:,d));
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
                 if size(find(abs(trainset_fitness(:,1) - archive_fitness(arch,1)) < 1e-3),1) == 0
                     trainset(train_num,:) = archive_position(arch,:);
                     trainset_fitness(train_num,1) = archive_fitness(arch,1);
                     train_num = train_num+1;
                 else
                     same = 1;
                     for t=1:train_num-1
                         a = find(abs(trainset(t,:)-archive_position(arch,:))>1e-1);
                         if size(a',1)>1
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
                 if size(find(abs(trainset_fitness(:,1) - archive_fitness(arch,1)) < 1e-3),1) == 0
                     trainset(train_num,:) = archive_position(arch,:);
                     trainset_fitness(train_num,1) = archive_fitness(arch,1);
                     train_num = train_num+1;
                 else
                     same = 1;
                     for t=1:train_num-1
                         a=find(abs(trainset(t,:)-archive_position(arch,:))>1e-1);
                         if size(a',1)>1
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
         
         if train_num>dimension*5
             [sortvalue sortid] = sort(trainset_fitness,'ascend');
             trainset1 = trainset;
             trainset_fitness1 = trainset_fitness;
             times = floor(train_num/dimension*5);
             selectedid = 1;
             for k=1:dimension*5
                 trainset(k,:) = trainset1(sortid(selectedid),:);
                 trainset_fitness(k,:) = trainset_fitness1(sortid(selectedid),:);
                 selectedid = selected+selectedid*times;
             end
         end
         
         for d=1:dimension
             max_range(d,1)=max(current_position(:,d));
             min_range(d,1)=min(current_position(:,d));
         end
         spread=0;
         for d=1:dimension
             spread=spread+abs(max_range(d,1)-min_range(d,1));
         end
         spread = sqrt(spread)*(mean(trainset_fitness))/dimension; 
         trainmodel = newrbe(trainset',trainset_fitness',spread);
         
         n_currarchive = n_archive-1;
         min_diff = abs(max(archive_fitness)-min(archive_fitness));
         max_diff = 0;
         mindiff_id = 0;
         maxdiff_id = 0;
         min_neighbordist = 0;
         min_neighborid = 0;
         min_distneighborvalue = min_diff;
         min_distneighborid = 0;
         pbest1=pbest;
         pbest_fitness1 = pbest_fitness;
         evaluated1 = evaluated;
         for i=1:popsize
             approximated_fitness(i,1) = sim(trainmodel,current_position(i,:)');
             current_fitness(i,1) = approximated_fitness(i,1);
             if current_fitness(i,:) < pbest_fitness(i,1)
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
                     min_neighborid = i;
                 end
             else
                 dist=zeros(n_archive-1,1);
                 for k=1:n_archive-1
                     temp_dist = 0;
                     for d=1:dimension
                         temp_dist = temp_dist + (current_position(i,d)-archive_position(k,d))^2;
                     end
                     dist(k,1) = sqrt(temp_dist);
                 end
                 [min_distvalue min_idvalue] = min(dist);
                 if min_distvalue>max_diff
                     max_diff = min_distvalue;
                     maxdiff_id = i;
                 end
             end
         end
         if min_neighborid ~= 0
             if func_id >= 9 && func_id <= 14
                 current_fitness(min_neighborid,1) = benchmark_func(current_position(min_neighborid,:),flag1);
             else
                 current_fitness(min_neighborid,1) = fitness(current_position(min_neighborid,:),func_id);
             end
             label_position(evaltimes+1,:) = current_position(min_neighborid,:);
             label_fitness(evaltimes+1,1) = current_fitness(min_neighborid,1);
             evaltimes=evaltimes+1;
             gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
             archive_position(n_archive,:) = current_position(min_neighborid,:);
             archive_fitness(n_archive,1) = current_fitness(min_neighborid,1);
             n_archive = n_archive + 1;
             if current_fitness(min_neighborid,1) < pbest_fitness(min_neighborid,1)
                 pbest(min_neighborid,:) = current_position(min_neighborid,:);
                 pbest_fitness(min_neighborid,1) = current_fitness(min_neighborid,1);
                 evaluated(min_neighborid,1) = 1;
             else
                 pbest(min_neighborid,:) = pbest1(min_neighborid,:);
                 pbest_fitness(min_neighborid,1) = pbest_fitness1(min_neighborid,1);
                 evaluated(min_neighborid,1) = evaluated1(min_neighborid,1);
             end
         else
             if maxdiff_id ~= 0
                 if func_id >= 9 && func_id <= 14
                     current_fitness(maxdiff_id,1) = benchmark_func(current_position(maxdiff_id,:),flag1);
                 else
                     current_fitness(maxdiff_id,1) = fitness(current_position(maxdiff_id,:),func_id);
                 end
                 label_position(evaltimes+1,:) = current_position(maxdiff_id,:);
                 label_fitness(evaltimes+1,1) = current_fitness(maxdiff_id,1);
                 evaltimes=evaltimes+1;
                 gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
                 archive_position(n_archive,:) = current_position(maxdiff_id,:);
                 archive_fitness(n_archive,1) = current_fitness(maxdiff_id,1);
                 n_archive = n_archive + 1;
                 if current_fitness(maxdiff_id,1) < pbest_fitness(maxdiff_id,1)
                     pbest(maxdiff_id,:) = current_position(maxdiff_id,:);
                     pbest_fitness(maxdiff_id,1) = current_fitness(maxdiff_id,1);
                     evaluated(maxdiff_id,1) = 1;
                 else
                     pbest(maxdiff_id,:) = pbest1(maxdiff_id,:);
                     pbest_fitness(maxdiff_id,1) = pbest_fitness1(maxdiff_id,1);
                     evaluated(maxdiff_id,1) = evaluated1(maxdiff_id,1);
                 end
             end
         end
         
         approx_diff = [];
         approx_diff(1:popsize,1) = abs(max(current_fitness));
         for i=1:popsize
             temp_trainset = trainset;
             temp_trainsetfit = trainset_fitness;
             if i~=min_neighborid
                 temp_trainset=[temp_trainset;current_position(i,:)];
                 temp_trainsetfit = [temp_trainsetfit; current_fitness(i,1)];
                 temp_trainmodel = newrbe(temp_trainset',temp_trainsetfit',spread);
                 position_fit = sim(temp_trainmodel,current_position(min_neighborid,:)');
                 approx_diff(i,1) = abs(position_fit - current_fitness(min_neighborid,1));
             else if i~=maxdiff_id
                     temp_trainset=[temp_trainset;current_position(i,:)];
                     temp_trainsetfit = [temp_trainsetfit; current_fitness(i,1)];
                     temp_trainmodel = newrbe(temp_trainset',temp_trainsetfit',spread);
                     position_fit = sim(temp_trainmodel,current_position(maxdiff_id,:)');
                     approx_diff(i,1) = abs(position_fit - current_fitness(maxdiff_id,1));
                 end
             end
         end
         temp_archive = [];
         for i=1:popsize
             if i~=min_neighborid 
                 if min_neighborid~=0
                     if (approx_diff(i,1) - abs(approximated_fitness(min_neighborid,1) - current_fitness(min_neighborid,1)))<-1e-3
                         archive_position(n_archive,:) = current_position(i,:);
                         archive_fitness(n_archive,1) = current_fitness(i,1);
                         n_archive = n_archive+1;
                         temp_archive = [temp_archive;i];
                     end
                 else
                     if i~=maxdiff_id 
                         if (approx_diff(i,1) - abs(approximated_fitness(maxdiff_id,1) - current_fitness(maxdiff_id,1)))<-1e-3
                             archive_position(n_archive,:) = current_position(i,:);
                             archive_fitness(n_archive,1) = current_fitness(i,1);
                             n_archive = n_archive+1;
                             temp_archive = [temp_archive;i];
                         end
                     end
                 end
             end
         end
         [minvalue minid] = min(pbest_fitness);
         if minvalue<gbest_fitness(1,1)
             foundid = find(temp_archive==minid);
             if size(foundid,1)>0
                 if minid == temp_archive(foundid,1)
                     if abs(current_fitness(minid,1)-pbest_fitness(minid,1))<1e-3 %& dist(current_position(minid,:),pbest(minid,:)')< 1e-3
                         if func_id >= 9 && func_id <= 14
                             pbest_fitness(minid,1) = benchmark_func(pbest(minid,:),flag1);
                         else
                             pbest_fitness(minid,1)=fitness(pbest(minid,:),func_id);
                         end
                         evaluated(minid,1) = 1;
                         label_position(evaltimes+1,:)=pbest(minid,:);
                         label_fitness(evaltimes+1,1)=pbest_fitness(minid,1);
                         evaltimes=evaltimes+1;
                         current_fitness(minid,1) = pbest_fitness(minid,1);
                         archive_position(n_archive-(size(temp_archive,1)-foundid+1),:) = pbest(minid,:);
                         archive_fitness(n_archive-(size(temp_archive,1)-foundid+1),1) = pbest_fitness(minid,:);
                         gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
                         if pbest_fitness(minid,1) < gbest_fitness
                             gbest(1,:) = pbest(minid,:);
                             gbest_fitness(1,1) = pbest_fitness(minid,1);
                         end
                     else
                         if evaluated(minid,1)==0
                             if func_id >= 9 && func_id <= 14
                                 pbest_fitness(minid,1) = benchmark_func(pbest(minid,:),flag1);
                             else
                                 pbest_fitness(minid,1)=fitness(pbest(minid,:),func_id);
                             end
                             evaluated(minid,1) = 1;
                             label_position(evaltimes+1,:)=pbest(minid,:);
                             label_fitness(evaltimes+1,1)=pbest_fitness(minid,1);
                             evaltimes=evaltimes+1;
                             archive_position(n_archive,:) = pbest(minid,:);
                             archive_fitness(n_archive,1) = pbest_fitness(minid,1);
                             n_archive = n_archive + 1;
                             gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
                             if pbest_fitness(minid,1) < gbest_fitness
                                 gbest(1,:) = pbest(minid,:);
                                 gbest_fitness(1,1) = pbest_fitness(minid,1);
                             end
                         else
                             gbest(1,:) = pbest(minid,:);
                             gbest_fitness(1,1) = pbest_fitness(minid,1);
                         end
                     end
                 end
             else
                 if evaluated(minid,1)==0
                     if func_id >= 9 && func_id <= 14
                         pbest_fitness(minid,1) = benchmark_func(pbest(minid,:),flag1);
                     else
                         pbest_fitness(minid,1)=fitness(pbest(minid,:),func_id);
                     end
                     label_position(evaltimes+1,:)=pbest(minid,:);
                     label_fitness(evaltimes+1,1)=pbest_fitness(minid,1);
                     evaluated(minid,1) = 1;
                     current_fitness(minid,1) = pbest_fitness(minid,1);
                     archive_position(n_archive,:) = pbest(minid,:);
                     archive_fitness(n_archive,1) = pbest_fitness(minid,1);
                     n_archive = n_archive + 1;
                     evaltimes=evaltimes+1;
                     gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
                     if pbest_fitness(minid,1) < gbest_fitness
                         gbest(1,:) = pbest(minid,:);
                         gbest_fitness(1,1) = pbest_fitness(minid,1);
                     end
                 else
                     gbest(1,:) = pbest(minid,:);
                     gbest_fitness(1,1) = pbest_fitness(minid,1);
                 end
             end
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