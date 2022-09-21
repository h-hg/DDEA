%% Core code of DESO
function [BP,BE,hisx,hisf,NFEs,CE,gfs,Strategy1,Strategy2]...
    = DESO(FUN, Dim, NP, L_bound, U_bound, initial_position,initial_fittness, hisx, hisf, NFEs, Max_NFEs, CE, sn1, gfs)
%% Parameters setting
%---- FUN       objective function
%---- Dim       Number of dimensions
%---- NP        Population Number
LB = repmat((L_bound),NP,1);              %---- Lower bound
UB = repmat((U_bound),NP,1);              %---- Upper bound
P = initial_position;                     % initialed Population position
E = initial_fittness;                     % initialed fittness
%---- hisx      history population position
%---- hisf      history sorted fittness
%---- NFEs      number of exact evaluation
%---- Max_NFEs  Maximum number of exact evaluation
%---- CE        achieve the exact fitness for plotting
%---- sn1       Compression parameter of convergence process record
%---- gfs       Sampling point according to fitness evaluation for ploting the convergence curve

BP = P(1,:);    % best population position
BE = E(1);      % best fittness

% DE parameter
F=0.5;
CR=0.9;

% number of local search samples
ls = 25 + Dim;
if ls > 60
    ls = 60;
end

% number of calling Strategys 
Strategy1 = 0;
Strategy2 = 0;

% show evaluated candidate each generation
show = 0; 

%% Main loop
while NFEs < Max_NFEs
    disp(['NFEs: ' num2str(NFEs)]); 
    
    Strategy = unidrnd(2);
    if Strategy == 1            % Strategy 1: Surrogate screening     
        Strategy1 = Strategy1 + 1; 
        LB = repmat((L_bound),length(E),1); 
        UB = repmat((U_bound),length(E),1);
        NP_last = 50; 
        NP_next = 50;
        P = hisx(1:NP_last,:); E = hisf(1:NP_last); % update P and E
        U = DEoperating(P,NP,Dim,hisx,F,CR,UB,LB);
        [ghf,id]=sort(hisf);                            % sort history data 
        gs=length(ghf(1:end));                          % global sample number  
        if gs > 300
            gs = 300;
        end
        ghx=hisx(id(1:gs),:);  ghf=ghf(1:gs);           % model_1 training samples
        ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
        spr=max(max(ghxd))/(Dim*gs)^(1/Dim);
        net=newrbe(ghx',ghf,spr);                       % newrbe - build RBF Neural Networks
        modelFUN=@(x) sim(net,x');                      % build global surrogate model
        fitnessModel=modelFUN(U);                       % surrogate model evaluate
        [fit,sidx]=sort(fitnessModel);                  % sort point based on fitness, get point indexs
        sam=U(sidx,:);                                  % sorted sample
        candidate_position = sam(1,:);                  % screening a candidate
    elseif  Strategy == 2       % Strategy 2: Surrogate sampling
        Strategy2 = Strategy2 + 1;
        [lhf, id]=sort(hisf);               
        lhx=hisx(id(1:ls),:); lhf = hisf(id(1:ls));
        lhxd=real(sqrt(lhx.^2*ones(size(lhx'))+ones(size(lhx))*(lhx').^2-2*lhx*(lhx')));
        spr=max(max(lhxd))/(Dim*ls)^(1/Dim);
        net=newrbe(lhx',lhf,spr);        % newrbe
        LocalModelFUN=@(x) sim(net,x');  % build local surrogate model 
        maxgen=100*Dim+1000; minerror=1e-20;
        [candidate_position,~] = DE(Dim, maxgen, LocalModelFUN, minerror, lhx); % find a optimum of surrogate model by optimizer
    end
    
    % judge Repeat Sample
    [~,ih,~]=intersect(hisx,candidate_position,'rows'); 
    if isempty(ih)~=1
        disp(['Sample Repeat and Delete it']);
        continue;
    end
    
    % evaluate candidate
    candidate_fit=FUN(candidate_position);
    NFEs = NFEs + 1; 
    if show 
        disp(['Strategy=' num2str(Strategy) ]);
        disp(['candidate_fit=' num2str(candidate_fit) ' point = ' num2str(candidate_position) ]);
    end
    
    % save candidate into dataset, and sort dataset
    hisx=[hisx; candidate_position];  hisf=[hisf, candidate_fit];   % update history database 
    [hisf,sidx]=sort(hisf);                                         % sort point based on fitness, get point indexs
    hisx=hisx(sidx,:); 
        
    % update BestEvaluation and BestPoint
    if  candidate_fit <= hisf(1) 
        BP = candidate_position;
        BE = candidate_fit;
        disp(['Best Cost(Strategy ' num2str(Strategy)  ') = ' num2str(BE) ' NFE=' num2str(NFEs)]);
    end
    
    % update CE for plotting
    CE(NFEs,:)=[NFEs,candidate_fit];
    if mod (NFEs,sn1)==0
        cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
    end
end