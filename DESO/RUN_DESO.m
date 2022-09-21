% DESO: Data-driven Evolutionary Sampling Optimization
function [ gsamp1 ,time_cost] = RUN_DESO(FUN, runs, Dim, L_bound, U_bound, fname) 
% 1.parameters
warning('off'); 
NP=50;                                      % Population size   
Max_NFEs=1000;                              % Maximum number of exact evaluation
CE=zeros(Max_NFEs,2);                       % Achive the exact fitness
sn1=1;                                      % Compression parameter of convergence process record
gfs=zeros(1,fix(Max_NFEs/sn1));             % Sampling point according to fitness evaluation for ploting the convergence curve
if Dim < 100                                % Initialization size
    initial_sample_size=100;     
elseif Dim >= 100                           
    initial_sample_size=200;
end

% 2.Runs
tic;
for r=1:runs
    % 2.1 Initialization 
    NFEs=0;                                         % counter of fitness evaluate
    fit=zeros(1,initial_sample_size);               % fitness
    sam=repmat(L_bound,initial_sample_size,1)+repmat((U_bound-L_bound),initial_sample_size,1).*lhsdesign(initial_sample_size,Dim);         % Initial LHS samples
    for i=1:initial_sample_size                                 
        fit(i)=FUN(sam(i,:));                                   
        NFEs=NFEs+1;
        CE(NFEs,:)=[NFEs,fit(i)];
        if mod(NFEs,sn1)==0
            cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
        end
    end
    hisx=sam; hisf=fit;                             % archive all the fitness evaluated samples        
    [r,NFEs]                                        % print  [r,NFEs]  
    
    % sort based on fitness
    [fit,sidx]=sort(fit);                           % sort point based on fitness, get point indexs
    sam=sam(sidx,:);                                % sorted sample
    initial_position=sam(1:NP,:);  initial_fittness=fit(1:NP);                % select best NP samples
    
    %--------------- Start search -------------------
    [pbest,pbestval,hisx,hisf,NFEs,CE,gfs,Strategy1,Strategy2]=DESO(FUN,Dim,NP,L_bound,U_bound,initial_position,initial_fittness,hisx,hisf,NFEs,Max_NFEs,CE,sn1,gfs);  
    %---------------- End search --------------------

    fprintf('Best fitness: %e\n',min(hisf));    
    gsamp1(r,:)=gfs;             					% record number of runs and optimization sample every time   
    run_time = toc;
    disp(run_time);
    usage(r,:)=[Strategy1,Strategy2];     
    cost_time(r,:)=run_time;     					% record time 
end    

% 3.output

% Evaluation index
best_samp   = min(gsamp1(:,end));
worst_samp  = max(gsamp1(:,end));
mean_samp   = mean(gsamp1(:,end));
median_samp = median(gsamp1(:,end));
std_samp    = std(gsamp1(:,end));
out1        = [best_samp,worst_samp,mean_samp,median_samp,std_samp];

% Convergence process (log10 and log are different!)
gsamp1_ave  = mean(gsamp1,1);
gsamp1_log  = log10(gsamp1_ave);
gsamplog    = log10(gsamp1);

% index Compressed convergence process
for j=1:Max_NFEs
    if mod(j,sn1)==0
        j1=j/sn1; gener_samp1(j1)=j;
    end
end
time_cost = cost_time
save(strcat('result/funcname/',num2str(Max_NFEs),'_',fname,' runs=',num2str(runs),' Dim=',num2str(Dim)));