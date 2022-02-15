clear all
cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));
run_times=1;
max_iteration=zeros(run_times,1);
evaluation_times=[];
evaluation_gbest=[];
hitcount_times=[];
best_fitness=[];
for i=1:run_times
    [evaluation_timestemp,evaluation_gbesttemp,max_itertemp,best_fit,hitcount_timestemp,global_approximate_output,local_approximate_output]=f_main_2_11_26;%f_main;
    for iter=1:max_itertemp
        evaluation_times(i,iter)=evaluation_timestemp(iter);
        evaluation_gbest(i,iter)=evaluation_gbesttemp(iter);
        approximate_global(i,iter)=global_approximate_output(iter);
        approximate_local(i,iter)=local_approximate_output(iter);
    end
    max_iteration(i,1)=max_itertemp;
    best_fitness(i,1)=best_fit;
end
save 'd:\f1_111.txt' evaluation_times -ascii;
save 'd:\f1_211.txt'  evaluation_gbest -ascii;
save 'd:\f1_311.txt' max_iteration -ascii;
save 'd:\f1_411.txt' best_fitness -ascii;
save 'd:\f1_511.txt' approximate_global -ascii;
save 'd:\f1_611.txt' approximate_local -ascii;

