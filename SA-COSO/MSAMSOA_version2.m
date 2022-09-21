function [gbest_output evaltimes_output final_best max_iteration N_evaluation G_optimum time_output] = MSAMOSA_version2(func_id)
global dimension
global popsize1
global popsize2
global max_length
global evaltimes

global initial_flag

initial_flag=0;
popsize1=30;
popsize2=200;
popsize = 230; %popsize = popsize1 + popsize2
run_times=1; %the times to run independently
max_length=1000; %maximum fitness evaluations
max_node=8; % number of nodes of the RBF hidden layer
c1=2.05;
phi=0.7298;
if (func_id==1 || func_id == 5)
    xmax=600; %Griewank
    xmin=-600;
    vmax=600;
    vmin=-600;
else if (func_id==2 || func_id== 6)
        xmax=32.768; %Ackley
        xmin=-32.768;
        vmax=32.768;
        vmin=-32.768;
    else if (func_id==3 || func_id == 7)
            xmax=2.048; %Rosenbrock
            xmin=-2.048;
            vmax=2.048;
            vmin=-2.048;
        else if (func_id==4 || func_id == 8)
                xmax=5.12; % ellipsoid
                xmin=-5.12;
                vmax=5.12;
                vmin=-5.12;
            else
                xmax=5;
                xmin=-5;
                vmax=5;
                vmin=-5;
            end
        end
    end
end
if (func_id <=4 || func_id == 9 || func_id == 11)
    dimension = 50;
else
    dimension = 100;
end

if (func_id == 9 || func_id == 10)
    global_flag = 10;
    initial_flag = 0;
else if (func_id == 11 || func_id == 12)
    global_flag = 19;
    initial_flag = 0;
    else
        global_flag = 0;
    end
end

max_archive=max_node * dimension + 10; %archive size = max_node * dimension + 10;
format short;
max_iteration=zeros(run_times,1);
gbest_output=[];
evaltimes_output=[];
time_output=[];

N_evaluation = [];
G_optimum = [];

for r=1:run_times
    initial_time=cputime;
    current_position=zeros(popsize1,dimension);
    current_velocity=zeros(popsize1,dimension);
    current_vel=zeros(popsize2,dimension);
    current_pos=zeros(popsize2,dimension);
    current_fit=zeros(popsize2,1);
    virtual_position=zeros(1,dimension);
    virtual_fitness=zeros(1,1);
    his_pos_1=zeros(popsize1,dimension);
    his_pos_2=zeros(popsize1,dimension);
    current_fitness=zeros(popsize1,1);
    his_fit_1=zeros(popsize1,1);
    his_fit_2=zeros(popsize1,1);
    pbest=zeros(popsize1,dimension);
    pbest_fitness=zeros(popsize1,1);
    rbfpbest=zeros(popsize2,dimension);
    rbfpbest_fitness=zeros(popsize2,1);
    dist_to_particle=zeros(popsize1,popsize1);
    r1=zeros(popsize1,dimension);
    r2=zeros(popsize1,dimension);
    r3=zeros(popsize1,dimension);
    fit_known=[];
    fit_eval=[];
    pbest_eval=zeros(popsize1,1);
    gbest=zeros(1,dimension);
    gbest_fitness=zeros(1,1);
    rbfgbest=zeros(1,dimension);
    rbfgbest_fitness=zeros(1,1);
    temp_save_pos=[];
    temp_save_fit=[];
    archive_pos=[];
    archive_fitness=[];
    rbfgbest_save=[];
    evaltimes=0;
    
    n_same = 0; %count the times that gbest is discarded
    
    initial_position = zeros(popsize,dimension); %generation of initial population including two swarms    
    %LHS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    latin_position = zeros(dimension, popsize);
    sample_selected = zeros(dimension, popsize);
    for i = 1:dimension
        for j = 1:popsize
            latin_position(i,j) = xmin + ((xmax-xmin)/popsize) * j; 
            sample_selected(i,j) = 0;
        end
    end
    
    pop_temp = popsize;
    
    for i = 1:popsize
        if(pop_temp == 1)
            for j = 1:dimension
                initial_position(i,j) = latin_position(j,1);
            end
        else
            latin_position_temp = zeros(dimension, pop_temp - 1);
            for j = 1:dimension
                selected_value = fix(rand() * pop_temp);
                while (selected_value == 0)
                    selected_value = fix(rand() * pop_temp);
                end
                initial_position(i,j) = latin_position(j,selected_value);
                a = latin_position(j,:);
                b = latin_position(j,selected_value);
                c = setdiff(a,b);
                latin_position_temp(j,:) = c;
            end
            latin_position = zeros(dimension, pop_temp);
            latin_position = latin_position_temp;
            pop_temp = pop_temp - 1;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate two populations
    for i=1:popsize1
        for t=1:dimension
            current_velocity(i,t)=rand*(vmax-vmin)+vmin;
            current_position(i,t)=initial_position(i,t);
        end
    end
    
    for i=1:popsize2
        for t=1:dimension
            current_vel(i,t)=rand*(vmax-vmin)+vmin;
            current_pos(i,t)=initial_position(popsize1+i,t);
        end
    end

    for i=1:popsize1
        temp_position=current_position(i,:);
        current_fitness(i,1)=function_run(temp_position,func_id,global_flag);
        evaltimes=evaltimes+1;
        fit_known(1,i,1)=1;
        fit_eval(1,i,1)=1;
    end
    for i=1:popsize2
        temp_position=current_pos(i,:);
        current_fit(i,1)=function_run(temp_position,func_id,global_flag);
        evaltimes=evaltimes+1;
    end
    
    %archive_pos and archive_fitness are used to save the information of
    %individuals that the fitness are evaluated using the real objective
    %function
    for i=1:popsize1
        archive_pos(i,:)=current_position(i,:);
        archive_fitness(i,1)=current_fitness(i,1);
    end
    for i=1:popsize2
        archive_pos(popsize1+i,:)=current_pos(i,:);
        archive_fitness(popsize1+i,1)=current_fit(i,1);
    end
    
    %define the value of parameters used in newrb, especially for spread,
    %as the value of spread will influence the performance of the RBF
    for t=1:dimension
        max_range(t,1)=max(archive_pos(:,t));
        min_range(t,1)=min(archive_pos(:,t));
    end
    spread=0;
    for t=1:dimension
        spread=spread+(max_range(t,1)-min_range(t,1))^2;
    end
    spread=sqrt(spread);
    spread_sum=spread;
    timenet=newrb(archive_pos',archive_fitness',0.1,spread,max_node,1);
    evaltimes_output(r,1)=evaltimes; % r represents the times of independencely run
    for i=1:popsize1
        pbest(i,:)=current_position(i,:);
        pbest_fitness(i,:)=current_fitness(i,:);
        pbest_eval(i,1)=1;
    end
    for i=1:popsize2
        rbfpbest(i,:)=current_pos(i,:);
        rbfpbest_fitness(i,1)=current_fit(i,1);
    end
    gbest(1,:)=pbest(1,:);
    gbest_fitness(1,1)=pbest_fitness(1,1);
    for i=2:popsize1
        if pbest_fitness(i,1)<gbest_fitness(1,1)
            gbest(1,:)=pbest(i,:);
            gbest_fitness(1,1)=pbest_fitness(i,1);
        end
    end
    rbfgbest(1,:)=rbfpbest(1,:);
    rbfgbest_fitness(1,1)=rbfpbest_fitness(1,1);
    
    for i=2:popsize2
        if rbfpbest_fitness(i,1)<rbfgbest_fitness(1,1)
            rbfgbest(1,:)=rbfpbest(i,:);
            rbfgbest_fitness(1,1)=rbfpbest_fitness(i,1);
        end
    end
    test_fit=min(gbest_fitness(1,1),rbfgbest_fitness(1,1))
        if test_fit==gbest_fitness(1,1)
            test_pos=gbest;
        else
            test_pos=rbfgbest;
        end
    gbest_output(r,1)=min(gbest_fitness(1,1),rbfgbest_fitness(1,1));
    totalupdate1=0;
    iter=1;
    rbfgbest_save(iter,1)=rbfgbest_fitness(1,1);
    
    N_evaluation(r,1) = popsize1+popsize2;
    G_optimum(r,1) = test_fit;
    evaltimes_old = evaltimes;
    
    while (evaltimes<1000) 
        save_id=0;
        c2 = 1.025;
        c3 = 1.025;
       
        for i=1:popsize2
            rbfuse_archive(i,:)=current_pos(i,:);
            rbfuse_archivefit(i,1)=current_fit(i,1);
        end
        
        selected_set = zeros(popsize2,1);
        for i=1:popsize2
            select_id=fix(rand*(size(archive_fitness,1)));
            while select_id==0
                select_id=fix(rand*(size(archive_fitness,1)));
            end
            chosen = 0;
            for j=1:size(selected_set,1)
                if (select_id == selected_set(j,1))
                    chosen = 1;
                    break;
                end  
            end
            while chosen == 1
                select_id=fix(rand*(size(archive_fitness,1)));
                while select_id==0
                    select_id=fix(rand*(size(archive_fitness,1)));
                end
                chosen = 0;
                for j=1:size(selected_set,1)
                    if (select_id == selected_set(j,1))
                        chosen = 1;
                        break;
                    end
                end
            end
            selected_set(i,1) = select_id;
            rbfuse_archive(popsize2+i,:)=archive_pos(select_id,:);
            rbfuse_archivefit(popsize2+i,1)=archive_fitness(select_id,1);
        end

        fit_known(iter+1,:,1)=0;
        fit_eval(iter+1,:,1)=0;
        % record the historical position and fitness of population 1;
        if iter==1
            his_pos_1=current_position;
            his_fit_1=current_fitness;
        else
            his_pos_2=his_pos_1;
            his_fit_2=his_fit_1;
            his_pos_1=current_position;
            his_fit_1=current_fitness;
        end
        [sorted_fit,sorted_index]=sort(his_fit_1);
        for i=1:popsize1 %update the velocity and position in population 1
            for t=1:dimension
                r1(i,t)=rand();
                r2(i,t)=rand();
                r3(i,t)=rand();
                current_velocity(i,t)=phi*(current_velocity(i,t)+c1*r1(i,t)*(pbest(i,t)-current_position(i,t))+c2*r2(i,t)*(gbest(1,t)-current_position(i,t))+c3*r3(i,t)*(rbfgbest(1,t)-current_position(i,t)));
                if current_velocity(i,t)>vmax
                    current_velocity(i,t)=vmax;
                end
                if current_velocity(i,t)<vmin
                    current_velocity(i,t)=vmin;
                end
                current_position(i,t)=current_position(i,t)+current_velocity(i,t);
                if current_position(i,t)>xmax
                    current_position(i,t)=xmax;
                end
                if current_position(i,t)<xmin
                    current_position(i,t)=xmin;
                end
            end
        end
        for i=1:popsize2
            %choose an individual to learn one dimension to follow
            better_ind=find(rbfuse_archivefit<current_fit(i,1));
            if size(better_ind,1)==1
                chose_id(1,:)=rbfuse_archive(better_ind(1,1),:);
            elseif size(better_ind,1)==2
                for t=1:dimension
                    if rand()>0.5
                        chose_id(1,t)=rbfuse_archive(better_ind(1,1),t);
                    else
                        chose_id(1,t)=rbfuse_archive(better_ind(2,1),t);
                    end
                end
            elseif size(better_ind,1)>2
                for t=1:dimension
                    generate_id=fix(rand*(size(better_ind,1)+1));
                    while generate_id==0 || generate_id==size(better_ind,1)+1
                        generate_id=fix(rand*(size(better_ind,1)+1));
                    end
                    chose_id(1,t)=rbfuse_archive(better_ind(generate_id,1),t);
                end
            end
            if size(better_ind,1)>0
                for t=1:dimension
                    current_vel(i,t)=rand*current_vel(i,t)+rand*(chose_id(1,t)-current_pos(i,t));
                    if current_vel(i,t)>vmax
                        current_vel(i,t)=vmax;
                    end
                    if current_vel(i,t)<vmin
                        current_vel(i,t)=vmin;
                    end
                    current_pos(i,t)=current_pos(i,t)+current_vel(i,t);
                    if current_pos(i,t)>xmax
                        current_pos(i,t)=xmax;
                    end
                    if current_pos(i,t)<xmin
                        current_pos(i,t)=xmin;
                    end
                end
                current_fit(i,1)=sim(timenet,current_pos(i,:)');
                rbfpbest(i,:)=current_pos(i,:);
                rbfpbest_fitness(i,1)=current_fit(i,1);
            end
        end
        fitness_determined = zeros(popsize1,1);
        fitness_determined(:,1) = 0;
        pop=1;
        hist_evaltimes = evaltimes;
        while pop<=popsize1
            i=sorted_index(pop,1);
            %i=pop;
            if iter==1
                temp_position=current_position(i,:);
                current_fitness(i,1)=function_run(temp_position,func_id,global_flag);
                temp_save_pos(save_id+1,:)=current_position(i,:);
                temp_save_fit(save_id+1,1)=current_fitness(i,1);
                save_id=save_id+1;
                evaltimes=evaltimes+1;
                N_evaluation(r,evaltimes-evaltimes_old+1) = evaltimes;
                G_optimum(r,evaltimes-evaltimes_old+1) = test_fit;
                fit_known(iter+1,i,1)=1;
                fit_eval(iter+1,i,1)=1;
            else
                if(fit_known(iter+1,i,1) == 0)
                    %use RBF to estimate
                    current_fitness(i,1) = sim(timenet,current_position(i,:)');
                    fit_known(iter+1,i,1) = 1;
                    fit_eval(iter+1,i,1) = 0;
                    fitness_determined(i,1) = 1;
                    for j=1:popsize1
                        if(i ~= j)
                            dist_to_particle(i,j) = distance_between(current_position(i,:),current_position(j,:),dimension);
                        else
                            dist_to_particle(i,i) = 100000; %set to the largest
                        end
                    end
                    ind=find(dist_to_particle(i,:)>0);
                    [min_dist,min_ind]=min(dist_to_particle(i,ind));
                    mi=sorted_index(ind(min_ind));
                    if((fit_known(iter+1,mi,1) == 0) || (fit_known(iter+1,mi,1) == 1 && fit_eval(iter+1,mi,1) == 0 && fitness_determined(mi,1) == 0))
                        for t=1:dimension
                            virtual_position(1,t)=current_position(i,t)+(1+phi-phi*c1*r1(mi,t)-phi*c2*r2(mi,t)-phi*c3*r3(mi,t))*his_pos_1(mi,t)+phi*his_pos_2(i,t)+phi*c1*r1(mi,t)*pbest(mi,t)+phi*c2*r2(mi,t)*gbest(1,t)+phi*c3*r3(mi,t)*rbfgbest(1,t);
                        end
                        dist_to_virtual(1,1) = distance_between(virtual_position(1,:),current_position(i,:),dimension);
                        dist_to_virtual(2,1) = distance_between(virtual_position(1,:),his_pos_1(mi,:),dimension);
                        dist_to_virtual(3,1) = distance_between(virtual_position(1,:),his_pos_2(i,:),dimension);
                        dist_to_virtual(4,1) = distance_between(virtual_position(1,:),pbest(mi,:),dimension);
                        dist_to_virtual(5,1) = distance_between(virtual_position(1,:),gbest(1,:),dimension);
                        dist_to_virtual(6,1) = distance_between(virtual_position(1,:),rbfgbest(1,:),dimension);
                        
                        dist_to_virtual(7,1) = distance_between(virtual_position(1,:),current_position(mi,:),dimension);
                        dist_to_virtual(8,1) = distance_between(virtual_position(1,:),his_pos_1(i,:),dimension);
                        dist_to_virtual(9,1) = distance_between(virtual_position(1,:),his_pos_2(mi,:),dimension);
                        dist_to_virtual(10,1) = distance_between(virtual_position(1,:),pbest(i,:),dimension);
                        dist_to_virtual(11,1) = dist_to_virtual(5,1);
                        dist_to_virtual(12,1) = dist_to_virtual(6,1);
                        if (size(find(dist_to_virtual == 0),1) == 0)
                            dist_temp1=1/dist_to_virtual(1,1)+1/dist_to_virtual(2,1)+1/dist_to_virtual(3,1)+1/dist_to_virtual(4,1)+1/dist_to_virtual(5,1)+1/dist_to_virtual(6,1);
                            dist_temp2=1/dist_to_virtual(7,1)+1/dist_to_virtual(8,1)+1/dist_to_virtual(9,1)+1/dist_to_virtual(10,1)+1/dist_to_virtual(11,1)+1/dist_to_virtual(12,1);
                            dist_radio = dist_temp2 / dist_temp1;
                            virtual_fitness(1,1) = current_fitness(i,1)/dist_to_virtual(1,1) + his_fit_1(mi,1)/dist_to_virtual(2,1) + his_fit_2(i,1)/dist_to_virtual(3,1) + pbest_fitness(mi,1)/dist_to_virtual(4,1) + gbest_fitness(1,1)/dist_to_virtual(5,1) + rbfgbest_fitness(1,1)/dist_to_virtual(6,1);
                            if(fit_known(iter+1,mi,1) == 1)
                                current_fitness(mi,1) = dist_to_virtual(7,1) * (virtual_fitness(1,1) * dist_radio - (his_fit_1(i,1)/dist_to_virtual(8,1)+his_fit_2(mi,1)/dist_to_virtual(9,1)+pbest_fitness(i,1)/dist_to_virtual(10,1)+gbest_fitness(1,1)/dist_to_virtual(11,1)+rbfgbest_fitness(1,1)/dist_to_virtual(12,1)));
                            else
                                temp_fit = dist_to_virtual(7,1) * (virtual_fitness(1,1) * dist_radio - (his_fit_1(i,1)/dist_to_virtual(8,1)+his_fit_2(mi,1)/dist_to_virtual(9,1)+pbest_fitness(i,1)/dist_to_virtual(10,1)+gbest_fitness(1,1)/dist_to_virtual(11,1)+rbfgbest_fitness(1,1)/dist_to_virtual(12,1)));
                                current_fitness(mi,1) = min(current_fitness(mi,1), temp_fit);
                            end
                            fit_known(iter+1,mi,1) = 1;
                            fit_eval(iter+1,mi,1) = 0;
                        end
                    end
                    ind=find(dist_to_particle(i,:)==0);
                    for k=1:size(ind',1)
                        if(fitness_determined(ind(1,k) == 0))
                            current_fitness(ind(1,k),1) = current_fitness(mi,1);
                        end
                    end
                else
                    fitness_determined(i,1) = 1;
                    for j=1:popsize1
                        if(i ~= j)
                            dist_to_particle(i,j) = distance_between(current_position(i,:),current_position(j,:),dimension);
                        else
                            dist_to_particle(i,i) = 100000; %set to the largest
                        end
                    end
                    ind=find(dist_to_particle(i,:)>0);
                    [min_dist,min_ind]=min(dist_to_particle(i,ind));
                    mi=sorted_index(ind(min_ind));
                    if((fit_known(iter+1,mi,1) == 0) || (fit_known(iter+1,mi,1) == 1 && fit_eval(iter+1,mi,1) == 0 && fitness_determined(mi,1) == 0))
                        for t=1:dimension
                            virtual_position(1,t)=current_position(i,t)+(1+phi-phi*c1*r1(mi,t)-phi*c2*r2(mi,t)-phi*c3*r3(mi,t))*his_pos_1(mi,t)+phi*his_pos_2(i,t)+phi*c1*r1(mi,t)*pbest(mi,t)+phi*c2*r2(mi,t)*gbest(1,t)+phi*c3*r3(mi,t)*rbfgbest(1,t);
                        end
                        dist_to_virtual(1,1) = distance_between(virtual_position(1,:),current_position(i,:),dimension);
                        dist_to_virtual(2,1) = distance_between(virtual_position(1,:),his_pos_1(mi,:),dimension);
                        dist_to_virtual(3,1) = distance_between(virtual_position(1,:),his_pos_2(i,:),dimension);
                        dist_to_virtual(4,1) = distance_between(virtual_position(1,:),pbest(mi,:),dimension);
                        dist_to_virtual(5,1) = distance_between(virtual_position(1,:),gbest(1,:),dimension);
                        dist_to_virtual(6,1) = distance_between(virtual_position(1,:),rbfgbest(1,:),dimension);
                        
                        dist_to_virtual(7,1) = distance_between(virtual_position(1,:),current_position(mi,:),dimension);
                        dist_to_virtual(8,1) = distance_between(virtual_position(1,:),his_pos_1(i,:),dimension);
                        dist_to_virtual(9,1) = distance_between(virtual_position(1,:),his_pos_2(mi,:),dimension);
                        dist_to_virtual(10,1) = distance_between(virtual_position(1,:),pbest(i,:),dimension);
                        dist_to_virtual(11,1) = dist_to_virtual(5,1);
                        dist_to_virtual(12,1) = dist_to_virtual(6,1);
                        if (size(find(dist_to_virtual == 0),1) == 0)
                            dist_temp1=1/dist_to_virtual(1,1)+1/dist_to_virtual(2,1)+1/dist_to_virtual(3,1)+1/dist_to_virtual(4,1)+1/dist_to_virtual(5,1)+1/dist_to_virtual(6,1);
                            dist_temp2=1/dist_to_virtual(7,1)+1/dist_to_virtual(8,1)+1/dist_to_virtual(9,1)+1/dist_to_virtual(10,1)+1/dist_to_virtual(11,1)+1/dist_to_virtual(12,1);
                            dist_radio = dist_temp2 / dist_temp1;
                            virtual_fitness(1,1) = current_fitness(i,1)/dist_to_virtual(1,1) + his_fit_1(mi,1)/dist_to_virtual(2,1) + his_fit_2(i,1)/dist_to_virtual(3,1) + pbest_fitness(mi,1)/dist_to_virtual(4,1) + gbest_fitness(1,1)/dist_to_virtual(5,1) + rbfgbest_fitness(1,1)/dist_to_virtual(6,1);
                            if(fit_known(iter+1,mi,1) == 1)
                                current_fitness(mi,1) = dist_to_virtual(7,1) * (virtual_fitness(1,1) * dist_radio - (his_fit_1(i,1)/dist_to_virtual(8,1)+his_fit_2(mi,1)/dist_to_virtual(9,1)+pbest_fitness(i,1)/dist_to_virtual(10,1)+gbest_fitness(1,1)/dist_to_virtual(11,1)+rbfgbest_fitness(1,1)/dist_to_virtual(12,1)));
                            else
                                temp_fit = dist_to_virtual(7,1) * (virtual_fitness(1,1) * dist_radio - (his_fit_1(i,1)/dist_to_virtual(8,1)+his_fit_2(mi,1)/dist_to_virtual(9,1)+pbest_fitness(i,1)/dist_to_virtual(10,1)+gbest_fitness(1,1)/dist_to_virtual(11,1)+rbfgbest_fitness(1,1)/dist_to_virtual(12,1)));
                                current_fitness(mi,1) = min(current_fitness(mi,1), temp_fit);
                            end
                            fit_known(iter+1,mi,1) = 1;
                            fit_eval(iter+1,mi,1) = 0;
                        end
                    end
                    ind=find(dist_to_particle(i,:)==0);
                    for k=1:size(ind',1)
                        if(fitness_determined(ind(1,k) == 0))
                            current_fitness(ind(1,k),1) = current_fitness(mi,1);
                        end
                    end
                end
            end
            pop=pop+1;
            if(current_fitness(i,1) ~=  sim(timenet,current_position(i,:)') && fit_eval(iter+1,i,1) == 0)
                if(current_fitness(i,1) < pbest_fitness(i,1) && sim(timenet,current_position(i,:)') < pbest_fitness(i,1))
                    temp_position=current_position(i,:);
                    current_fitness(i,1)=function_run(temp_position,func_id,global_flag);
                    temp_save_pos(save_id+1,:)=current_position(i,:);
                    temp_save_fit(save_id+1,1)=current_fitness(i,1);
                    save_id=save_id+1;
                    evaltimes=evaltimes+1;
                    N_evaluation(r,evaltimes-evaltimes_old+1) = evaltimes;
                    G_optimum(r,evaltimes-evaltimes_old+1) = test_fit;
                    fit_known(iter+1,i,1)=1;
                    fit_eval(iter+1,i,1)=1;
                end
            end
        end
        
        %using the second method
        if(hist_evaltimes == evaltimes)
            %%%calculate the mean error
            avg_err = 0;
            num = 0;
            for i=1:popsize1
                if current_fitness(i,1) ~= sim(timenet,current_position(i,:)')
                    avg_err = avg_err + abs(current_fitness(i,1) - sim(timenet,current_position(i,:)'));
                    num = num + 1;
                end
            end
            avg_err = avg_err / num;
            for i=1:popsize1
                if current_fitness(i,1) ~= sim(timenet,current_position(i,:)')
                    if abs(current_fitness(i,1) - sim(timenet,current_position(i,:)')) > avg_err
                        temp_position=current_position(i,:);
                        current_fitness(i,1)=function_run(temp_position,func_id,global_flag);
                        temp_save_pos(save_id+1,:)=current_position(i,:);
                        temp_save_fit(save_id+1,1)=current_fitness(i,1);
                        save_id=save_id+1;
                        evaltimes=evaltimes+1;
                        N_evaluation(r,evaltimes-evaltimes_old+1) = evaltimes;
                        G_optimum(r,evaltimes-evaltimes_old+1) = test_fit;
                        fit_known(iter+1,i,1)=1;
                        fit_eval(iter+1,i,1)=1;
                    end
                end
            end
        end
        
        for i=1:popsize1
            if current_fitness(i,1)<pbest_fitness(i,1)
                pbest(i,:)=current_position(i,:);
                pbest_fitness(i,1)=current_fitness(i,1);
                pbest_eval(i,1)=fit_eval(iter+1,i,1);
            end
        end
        
        best_temp = pbest_fitness(1,1);
        best_id = 1;
        for i=2:popsize1
            if pbest_fitness(i,1) < best_temp
                best_temp = pbest_fitness(i,1);
                best_id = i;
            end
        end
        if pbest_eval(best_id,1)== 1
            if(gbest_fitness(1,1) > pbest_fitness(best_id,1))
                gbest(1,:)=pbest(best_id,:);
                gbest_fitness(1,1)=pbest_fitness(best_id,1);
            end
        else
            temp_position=pbest(best_id,:);
            pbest_fitness(best_id,1)=function_run(temp_position,func_id,global_flag);
            evaltimes=evaltimes+1;
            N_evaluation(r,evaltimes-evaltimes_old+1) = evaltimes;
            G_optimum(r,evaltimes-evaltimes_old+1) = test_fit;
            temp_save_pos(save_id+1,:)=pbest(best_id,:);
            temp_save_fit(save_id+1,1)=pbest_fitness(best_id,1);
            save_id=save_id+1;
            pbest_eval(best_id,1)=1;
            if pbest_fitness(best_id,1)<gbest_fitness(1,1)
                gbest(1,:)=pbest(best_id,:);
                gbest_fitness(1,1)=pbest_fitness(best_id,1);
            end
        end
        
        rbfbest_temp = rbfpbest_fitness(1,1);
        rbfbest_id = 1;
        for i=2:popsize2
            if rbfpbest_fitness(i,1) < rbfbest_temp
                rbfbest_temp = rbfpbest_fitness(i,1);
                rbfbest_id = i;
            end
        end
        temp_position=rbfpbest(rbfbest_id,:);
        rbfpbest_fitness(rbfbest_id,1)=function_run(temp_position,func_id,global_flag);
        current_fit(rbfbest_id,1)=rbfpbest_fitness(rbfbest_id,1);
        temp_save_pos(save_id+1,:)=rbfpbest(rbfbest_id,:);
        temp_save_fit(save_id+1,1)=rbfpbest_fitness(rbfbest_id,1);
        save_id=save_id+1;
        evaltimes=evaltimes+1;
        N_evaluation(r,evaltimes-evaltimes_old+1) = evaltimes;
        G_optimum(r,evaltimes-evaltimes_old+1) = test_fit;
        if rbfpbest_fitness(rbfbest_id,1)<rbfgbest_fitness(1,1)
            rbfgbest(1,:)=rbfpbest(rbfbest_id,:);
            rbfgbest_fitness(1,1)=rbfpbest_fitness(rbfbest_id,1);
        end
       
        gbest_output(r,iter+1)=min(gbest_fitness(1,1),rbfgbest_fitness(1,1));
        [gbest_fitness(1,1) rbfgbest_fitness(1,1)]
        test_fit=min(gbest_fitness(1,1),rbfgbest_fitness(1,1))
        if test_fit==gbest_fitness(1,1)
            test_pos=gbest;
        else
            test_pos=rbfgbest;
        end
        evaltimes_output(r,iter+1)=evaltimes;
        iter=iter+1;
        evaltimes
        rbfgbest_save(iter,1)=rbfgbest_fitness(1,1);
        %update the RBF network
        
        dist_to_each=[];
        dist_to_archive=[];
        min_each=[];
        min_each_id=[];
        max_each=[];
        max_each_id=[];
        dist_to_pop=[];

          %using the mean error of the individuals in the temp_save to
          %decide whether the individual should be saved or not
          meanfit_eval_esti = 0;
          for k1=1:size(temp_save_fit,1)
              meanfit_eval_esti = meanfit_eval_esti + abs(temp_save_fit(k1,1) - sim(timenet,temp_save_pos(k1,:)'));
          end
          meanfit_eval_esti = meanfit_eval_esti/size(temp_save_fit,1);
          for k1=1:size(archive_fitness,1) %calculate the distance between each individual in the archive and the individual in the population2
            for k2=1:size(current_fit,1)
                archive_to_pop(k1,k2) = distance_between(archive_pos(k1,:),current_pos(k2,:),dimension);
            end
            minarc_to_pop(k1,1)=min(archive_to_pop(k1,:));
        end
        [max_minarc_pop max_minarc_popid]=max(minarc_to_pop);
        %final revision
        %add judgement whether the gbest position of both population will
        %be discarded
        
        dist_archive = zeros(size(temp_save_fit,1),size(archive_fitness,1));
        for k1=1:size(temp_save_fit,1)
            for k2=1:size(archive_fitness,1)
                dist_archive(k1,k2) = distance_between(temp_save_pos(k1,:),archive_pos(k2,:),dimension);
            end
            same=find(dist_archive(k1,:) < 10^(-4));
            if size(same',1)==0
                    if size(archive_fitness,1)<max_archive
                        archive_pos(size(archive_pos,1)+1,:)=temp_save_pos(k1,:);
                        archive_fitness(size(archive_fitness,1)+1,1)=temp_save_fit(k1,1);
                    else
                        for k2=1:size(current_fit,1)
                            dist_to_pop(k1,k2) = distance_between(temp_save_pos(k1,:),current_pos(k2,:),dimension);
                        end
                        min_to_pop(k1,1)=min(dist_to_pop(k1,:));
                        if(min_to_pop(k1,1) < max_minarc_pop)
                            same_to_gbest = 1;
                            for d=1:dimension
                                if abs(archive_pos(max_minarc_popid,d) - gbest(1,d)) > 1e-7
                                    same_to_gbest = 0;
                                    break;
                                end
                            end
                            if same_to_gbest == 1
                                n_same = n_same + 1;
                            end
                            
                            archive_pos(max_minarc_popid,:)=temp_save_pos(k1,:);
                            archive_fitness(max_minarc_popid,1)=temp_save_fit(k1,1);
                            archive_to_pop(max_minarc_popid,:)=dist_to_pop(k1,:);
                            minarc_to_pop(max_minarc_popid,1)=min_to_pop(k1,1);
                            [max_minarc_pop max_minarc_popid]=max(minarc_to_pop);
                        else
                            %%%if the distance is larger
                            same_to_gbest = 1;
                            for d=1:dimension
                                if abs(temp_save_pos(k1,d)-gbest(1,d)) > 1e-7
                                    same_to_gbest = 0;
                                    break;
                                end
                            end
                            if same_to_gbest == 1
                                n_same = n_same + 1;
                            end
                        end
                    end
                %end
            end
        end
        
        for t=1:dimension
            max_range(t,1)=max(archive_pos(:,t));
            min_range(t,1)=min(archive_pos(:,t));
        end

        spread1=0;
        for t=1:dimension
            spread1=spread1+(max_range(t,1)-min_range(t,1))^2;
        end
        spread1=sqrt(spread1);
         %spread=spread1/max_node
        spread_sum = spread_sum+spread1;
        spread = spread_sum/(iter+1);
        
        timenet=newrb(archive_pos',archive_fitness',0.1,spread,max_node,1);

    end
    N_times(r,1) = n_same;
    time_output(r,1) = cputime - initial_time;
    final_best(r,1)=min(gbest_fitness(1,1),rbfgbest_fitness(1,1));
    max_iteration(r,1)=iter-1
    max_node
    
end

% if func_id==1
%     %griewank\
% save 'd:\collision1.txt' N_times -ascii;
% else if func_id==2
%         %\ackley
%         save 'd:\collision2.txt' N_times -ascii;
%     else if func_id==3
%             %\rosenbrock
%             save 'd:\collision3.txt' N_times -ascii;
%         else if func_id==4
%                 %\\ellipsoid
%                 save 'd:\collision4.txt' N_times -ascii;
%             end
%         end
%     end
% end

final_best
min_value=min(final_best)
mean_value=mean(final_best)
max_value=max(final_best)
end