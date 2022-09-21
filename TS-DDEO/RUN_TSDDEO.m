%------------------------------------------------------------------------
% This code is part of the program that produces the results in the following paper:
% Huixiang Zhen, Wenyin Gong, Ling Wang, Fei Ming, and Zuowen Liao. "Two-stage Data-driven Evolutionary Optimization for High-dimensional Expensive Problems", IEEE Transactions on Cybernetics, accepted, 2021.
% You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%----------------------------------------------------------------------------------------------------------------------------------------

function [ gsamp1 ,time_cost] = RUN_TSDDEO(runs, D, FUN, LB, UB, fname)
time_begin=tic;
warning('off');
addpath(genpath(pwd));

for r=1:runs
    % main loop
    fprintf('\n');
    disp(['FUNCTION: ', fname,' RUN: ', num2str(r)]);  
    fprintf('\n');
    [hisx,hisf,fitcount,mf,CE,gfs]= TSDDEO(FUN,D,LB,UB); 
    fprintf('Best fitness (PSO-final): %e\n',min(hisf));       
    gsamp1(r,:)=gfs(1:mf);
end    

%%%%%%%%%%%%%%%%%%%%% Output options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% usage_ave=mean(usage); 
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
time_cost=toc(time_begin);
save(strcat('result/NFE',num2str(mf),'_',fname,' runs=',num2str(runs),' Dim=',num2str(D)));
end
