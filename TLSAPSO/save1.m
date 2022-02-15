save 'd:\a.txt' evaltimes_output -ascii;
save 'd:\b.txt' evaltimes_output1 -ascii;
save 'd:\c.txt' evaluation_gbest -ascii;

save 'd:\a.txt' evaluation_times -ascii;
save 'd:\c.txt'  evaluation_gbest -ascii;

save 'd:\f1_times.txt' evaluation_times -ascii;
save 'd:\f1_gbest.txt'  evaluation_gbest -ascii;
a=zeros(5,1);
for i=1:run_times
    a(i,1)=evaluation_gbest(i,max_iteration(i));
end
save 'd:\gbest_f1.txt' a -ascii;
save 'd:\iteration_f1.txt' max_iteration -ascii;

save 'd:\f22_times.txt' evaluation_times -ascii;
save 'd:\f22_gbest.txt'  evaluation_gbest -ascii;
a=zeros(5,1);
for i=1:run_times
    a(i,1)=evaluation_gbest(i,max_iteration(i));
end
save 'd:\gbest_f22.txt' a -ascii;
save 'd:\iteration_f22.txt' max_iteration -ascii;

save 'd:\f22_times.txt' evaluation_times -ascii;
save 'd:\f22_gbest.txt'  evaluation_gbest -ascii;
a=zeros(5,1);
for i=1:run_times
    a(i,1)=evaluation_gbest(i,max_iteration(i));
end
save 'd:\gbest_f22.txt' a -ascii;
save 'd:\iteration_f22.txt' max_iteration -ascii;

save 'd:\dpso_10_2_times.txt' evaluation_times -ascii;
save 'd:\dpso_10_2_gbest.txt'  evaluation_gbest -ascii;
a=zeros(5,1);
for i=1:run_times
    a(i,1)=evaluation_gbest(i,max_iteration(i));
end
save 'd:\dpso_10_2_gbestfit.txt' a -ascii;
save 'd:\dpso_10_2_iteration.txt' max_iteration -ascii;