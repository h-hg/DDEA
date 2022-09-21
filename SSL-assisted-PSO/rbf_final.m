function [G_Optimal Iter_OptimalResult Time] = SPSO_w(func_id,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations)
%function [G_Optimal Time] = SPSO_w(func_id,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations)
% clc
% clear all
% global dimension
% global popsize
warning off
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
maxnode=8;
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
    labeled_position = [];
    labeled_fitness = [];
    labeled_iter = [];
    archive_iter = [];
    n_label = 1;
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
        labeled_position(n_label,:) = current_position(i,:);
        labeled_fitness(n_label,1) = current_fitness(i,1);
        labeled_iter(n_label,1) = iter;
        archive_iter(n_archive,1) = iter;
        n_label = n_label + 1;
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
         iter = iter +1;
         [sort_pbestfit sort_pbestid] = sort(pbest_fitness);
         avg = mean(current_position);
         for i=1:popsize
             sortid = find(sort_pbestid == i);
             for d=1:dimension
                 randset = randperm(sortid);
                 selected_pbest = sort_pbestid(randset(1,1));
                 current_velocity(i,d) = rand*current_velocity(i,d) + rand*(pbest(selected_pbest,d)-current_position(i,d))+0.001*rand*(avg(1,d)-current_position(i,d));
%                  current_velocity(i,d) = 0.7293*(current_velocity(i,d) + 2.05*rand()*(pbest(i,d) - current_position(i,d)) + 2.05*rand()*(gbest(1,d) - current_position(i,d)));
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
         end %update the velocity and position
         
         setsize = 2*(dimension+1);%

         labeltrain_num = 1;
         labeltrainset = [];
         labeltrainset_fitness = [];
         [labelsortvalue labelsortid] = sort(labeled_fitness,'ascend');%sort(labeled_iter,'descend');%
         for arch1=1:n_label-1
             if labeltrain_num-1<setsize
                 arch = labelsortid(arch1);
                 if size(labeltrainset_fitness,1)==0
                     labeltrainset(labeltrain_num,:) = labeled_position(arch,:);
                     labeltrainset_fitness(labeltrain_num,1) = labeled_fitness(arch,1);
                     labeltrain_num = labeltrain_num+1;
                 else
                     if size(find(abs(labeltrainset_fitness(:,1) - labeled_fitness(arch,1)) < 1e-3),1) == 0
                         labeltrainset(labeltrain_num,:) = labeled_position(arch,:);
                         labeltrainset_fitness(labeltrain_num,1) = labeled_fitness(arch,1);
                         labeltrain_num = labeltrain_num+1;
                     else
                         same = 1;
                         for t=1:labeltrain_num-1
                             a=find(abs(labeltrainset(t,:)-labeled_position(arch,:))>(xmax-xmin)*0.01);
                             if size(a',1)>1
                                 same = 0;
                                 break;
                             end
                         end
                         if same == 0
                             labeltrainset(labeltrain_num,:) = labeled_position(arch,:);
                             labeltrainset_fitness(labeltrain_num,1) = labeled_fitness(arch,1);
                             labeltrain_num = labeltrain_num+1;
                         end
                     end
                 end
             end
         end
         labeltrain_num = labeltrain_num-1;
         
         flag ='cubic';
        [labellambda,labelgamma] = RBF(labeltrainset,labeltrainset_fitness,flag);
        label_FUN = @(x)RBF_eval(x,labeltrainset,labellambda,labelgamma,flag);
        
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
         newavg = mean(current_position);
         dist = zeros(popsize,1);
         for i=1:popsize
             dist(i,1) = 0;
             for d=1:dimension
                 dist(i,1) = dist(i,1)+(current_position(i,d)-newavg(1,d))^2;
             end
             dist(i,1) = sqrt(dist(i,1));
         end
         [centervalue centerid] = min(dist);
         
         labelapproximated_fitness = label_FUN(current_position);
         current_approximated = ones(popsize,1);      
         hasevaluated = 0;
         for i=1:popsize
             if labelapproximated_fitness(i,1)<pbest_fitness
                 if func_id >= 9 && func_id <= 14
                     current_fitness(i,1) = benchmark_func(current_position(i,:),flag1);
                 else
                     current_fitness(i,1) = fitness(current_position(i,:),func_id);
                 end
                 hasevaluated = 1;
                 current_approximated(i,1) = 0;
                 evaltimes=evaltimes+1;
                 gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
                 labeled_position(n_label,:) = current_position(i,:);
                 labeled_fitness(n_label,1) = current_fitness(i,1);
                 labeled_iter(n_label,1) = iter;
                 n_label = n_label+1;
                 if current_fitness(i,1)<pbest_fitness(i,1)
                     pbest(i,:) = current_position(i,:);
                     pbest_fitness(i,1) = current_fitness(i,1);
                     evaluated(i,1) = 1;
                 end
             end
         end
         if hasevaluated == 0
             [mincurrentvalue mincurrentid] = min(current_fitness);
             if func_id >= 9 && func_id <= 14
                 current_fitness(i,1) = benchmark_func(current_position(i,:),flag1);
             else
                 current_fitness(i,1) = fitness(current_position(i,:),func_id);
             end
             hasevaluated = 1;
             current_approximated(i,1) = 0;
             evaltimes=evaltimes+1;
             gbest_output(evaltimes-popsize+1) = gbest_fitness(1,1);
             labeled_position(n_label,:) = current_position(i,:);
             labeled_fitness(n_label,1) = current_fitness(i,1);
             labeled_iter(n_label,1) = iter;
             n_label = n_label+1;
             if current_fitness(i,1)<pbest_fitness(i,1)
                 pbest(i,:) = current_position(i,:);
                 pbest_fitness(i,1) = current_fitness(i,1);
                 evaluated(i,1) = 1;
             end
         end
         
         [minvalue minid] = min(pbest_fitness);
         if minvalue<gbest_fitness(1,1)
             if evaluated(minid,1) == 0
                 if func_id >= 9 && func_id <= 14
                     pbest_fitness(minid,1) = benchmark_func(pbest(minid,:),flag1);
                 else
                     pbest_fitness(minid,1)=fitness(pbest(minid,:),func_id);
                 end
                 evaluated(minid,1) = 1;
                 current_fitness(minid,1) = pbest_fitness(minid,1);
                 evaltimes=evaltimes+1;
                 labeled_position(n_label,:) = pbest(minid,:);
                 labeled_fitness(n_label,1) = pbest_fitness(minid,1);
                 labeled_iter(n_label,1) = iter;
                 n_label = n_label + 1;
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
         %iter = iter + 1
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

%use the center individual to judge whether it should be saved to archive