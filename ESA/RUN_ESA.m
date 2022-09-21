% Run ESA
function [ gsamp1 ,time_cost] = RUN_ESA(runs, D, FUN, LB, UB, fname)
time_begin = tic;
for r = 1:runs
    % main loop
    fprintf('\n');
    disp(['Fname:', fname,'  Run:', num2str(r)]);  
    fprintf('\n');
    [hisx,hisf,fitcount,mf,CE,gfs] = ESA(FUN,D,LB,UB); 
    fprintf('Best fitness: %e\n',min(hisf));
    gsamp1(r,:) = gfs(1:mf);
end    

%%%%%%%%%%%%%%%%%%%%% Output options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
best_samp   = min(gsamp1(:,end));
worst_samp  = max(gsamp1(:,end));
samp_mean   = mean(gsamp1(:,end));
samp_median = median(gsamp1(:,end));
std_samp    = std(gsamp1(:,end));
out1        = [best_samp,worst_samp,samp_mean,samp_median,std_samp];
gsamp1_ave  = mean(gsamp1,1);
gsamp1_log  = log10(gsamp1_ave);   
gsamplog    = log10(gsamp1);  

% Time Complexity
time_cost = toc(time_begin);
time_cost = time_cost/runs;
save(strcat('result/NFE',num2str(mf),'_',fname,' runs=',num2str(runs),' Dim=',num2str(D)));
end
