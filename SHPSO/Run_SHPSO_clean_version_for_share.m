%%%*********************************************************************************************%%%
%% Implementation of Surrogate-assisted Hierarchical Particle Swarm Optimization (SHPSO)
%% H. Yu, Y. Tan, J. Zeng, C. Sun, Y. Jin, Surrogate-assisted hierarchical 
%% particle swarm optimization, Information Sciences, 454-455 (2018) 59-72.
%%%*********************************************************************************************%%%
%% This paper and this code should be referenced whenever they are used to 
%% generate results for the user's own research. 
%%%*********************************************************************************************%%%
%% This matlab code was written by Haibo Yu
%% Please refer with all questions, comments, bug reports, etc. to tyustyuhaibo@126.com
% 
clear,clc,
time_begin=tic;
% rand('seed',sum(100*clock));
warning('off');

variable_domain; % Test functions' variable domain

D=30;       % dimension -----> user defined: 30 / 50 / 100 /
mf=1000;    % maximum number of exact evaluation -----> user defined: 1000(30D,50D,100D) / 8000(30D) /
% 
sn1=2;  
gfs=zeros(1,fix(mf/sn1));   %-------sampling point according to fitness evaluation for ploting the convergence curve
% 
pop_size=50;                % population size                
CE=zeros(mf,2);             % achive the exact fitness
% 
%% User defined function category
% fname='FITNESS';                      % Basic Function
fname='benchmark_func';             % 2005 CEC Benchmark Function: P. N. Suganthan et al. Problem Definitions and Evaluation Criteria for the CEC 2005 Special Session on Real-Parameter Optimization
% 
if D < 100
    initial_sample_size=100;     % < 100 dimension
elseif D == 100
    initial_sample_size=200;      % = 100 dimension
end
% 
runs=20;
for r=1:runs
    % Initialization procedure
    fitcount=0;
    sam=Xmin+(Xmax-Xmin).*lhsdesign(initial_sample_size,D);      % initial LHS samples
    fit=zeros(1,initial_sample_size);
    for i=1:initial_sample_size
        fit(i)=feval(fname,sam(i,:));                            % fitness of initial sample set
        fitcount=fitcount+1;
        CE(fitcount,:)=[fitcount,fit(i)];
        if mod (fitcount,sn1)==0
            cs1=fitcount/sn1; gfs(1,cs1)=min(CE(1:fitcount,2));
        end
    end
    hisx=sam; hisf=fit;                                          % archive all the exact evaluated samples
%    
    % main loop
    [r,fitcount]         
    [~,sidx]=sort(fit);                                         % 对个体按适应值由小到大进行排序
    sam=sam(sidx,:);  fit=fit(sidx);                            % 排序后的样本
    lox=sam(1:pop_size,:);  lof=fit(1:pop_size);                % 选择前ss个最优样本
    % -------------- P S O search --------------------- 
    psam=lox; efit=lof;                                         % PSO初始种群
    LB=repmat(Xmin,1,D);UB=repmat(Xmax,1,D);                    % 全空间上下界
    ps=size(psam,1);                                            % 种群规模
    [fpso,hisx,hisf,fitcount,CE,gfs,pbest,pbestval,formula1,formula2]= SHPSO_clean_version_for_share(fname,D,ps,LB,UB,psam,efit,hisx,hisf,fitcount,mf,CE,sn1,gfs); 
    %------------End P S O search-----------------------
% 
    fprintf('Best fitness (PSO-final): %e\n',min(hisf));        % 注：这个是采样点，并非当前最优点
% 
    gsamp1(r,:)=gfs;
    
%     usage(r,:)=[formula1,formula2]; 

end    
% usage_ave=mean(usage); % for SHPSO
%%%%%%%%%%%%%%%%%%%%% Output options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
best_samp   = min(gsamp1(:,end));
worst_samp  = max(gsamp1(:,end));
samp_mean   = mean(gsamp1(:,end));
samp_median = median(gsamp1(:,end));
std_samp    = std(gsamp1(:,end));
out1        = [best_samp,worst_samp,samp_mean,samp_median,std_samp];
gsamp1_ave  = mean(gsamp1,1);
gsamp1_log  = log(gsamp1_ave);
gsamplog    = log(gsamp1);  % for ackley function
for j=1:mf
    if mod(j,sn1)==0
        j1=j/sn1; gener_samp1(j1)=j;
    end
end
% 
figure(1);
plot(gener_samp1,gsamp1_log,'.-k','Markersize',16)
legend('SHPSO'); % xlim([50,1000]);
xlabel('Function Evaluation Calls ');
ylabel('Mean Fitenss Value (Natural log)');
% title('Shifted Rotated Rastrigin Function (F10)')   % cec05func * 10 * for 30 & 50 & 100 dimension
% title('Shifted Rotated Rastrigin Function (F16)')   % cec05func * 16 * for 30 & 50 & 100 dimension
% title('Rotated Hybrid Composition Function (F19)')  % cec05func * 19 * for 30 & 50 & 100 dimension
% title('Ackley Function')
title('Griewank Function')
% title('Rastrigin Function')
% title('Rosenbrock Function')
% title('Ellipsoid Function')

% Time Complexity
time_cost=toc(time_begin);

