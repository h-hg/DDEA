% This is code of ESA written by Huixiang Zhen, please refer all questions, comments, bug reports, etc. to zhenhuixiang@cug.edu.cn
function [hx,hf,NFE,MaxNFE,Archive_FEs,Archive_convergence] = ESA(FUN,Dim,L_bound,U_bound)
% Parameters setting
NFE = 0;
MaxNFE = 1000;                    
Archive_FEs = zeros(MaxNFE,2);  
Archive_convergence = zeros(1,MaxNFE);   
show = 0; % show evaluated candidate each FE

% Initial LHS samples
if Dim < 100
    initial_sample_size = 100;    
elseif Dim >= 100
    initial_sample_size = 150;                
end
sam = repmat(L_bound,initial_sample_size,1)+(repmat(U_bound,initial_sample_size,1)-repmat(L_bound,initial_sample_size,1)).*lhsdesign(initial_sample_size,Dim);
fit = zeros(1,initial_sample_size);
for i=1:initial_sample_size
    fit(i) = FUN(sam(i,:)); 
    NFE = NFE+1;
    Archive_FEs(NFE,:) = [NFE,fit(i)]; 
    Archive_convergence(1,i) = min(Archive_FEs(1:NFE,2));
end

% Build database
hx = sam; hf = fit;                                             
[~,sidx] = sort(hf);                                         
hx = hx(sidx,:);  hf = fit(sidx);  % history data

% DE parameters for action 1
NP = 50;
F=0.5;
CR=0.9;

% Samples number of local surrogate models for action 2-4
ls = 25 + Dim;  ls = min([ls,60]);
ls2 = 100; 

% RL parameters and initialization
alp = 0.1;                   
gamma = 0.9;   
Q_Agent = [ 0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;];
State = 1; 
Action = 1; 
        
% number of calling sampling Actions 
Action1 = 0;
Action2 = 0;
Action3 = 0;
Action4 = 0;
        
% Main loop
while NFE <= MaxNFE
    % Update state and action 
    R = 0; action_success = 0;
    Qvalue1 = Q_Agent(State,:); 
    temp = exp(Qvalue1);
    ratio = cumsum(temp)/sum(temp); 
    jtemp = find(rand(1)<ratio);
    Action = jtemp(1); 
    
    % Agent log 
    log_Q_Agent{NFE} = Q_Agent;
    log_State{NFE} = State;
    log_ratio{NFE} = ratio;
    log_Action{NFE} = Action;
    disp(['NFE: ' num2str(NFE) ' Action:' num2str(Action)] ); 
    
    % Execute action and obtain new data
    if Action == 1         
        % Action 1: Surrogate screening     
        Action1 = Action1 + 1; 
        LB = repmat((L_bound),NP,1); 
        UB = repmat((U_bound),NP,1);
        P = hx(1:NP,:); 
        U = DEoperating(P,NP,Dim,hx,F,CR,UB,LB);
        % build global surrogate model
        [ghf,id] = sort(hf);                            
        gs = length(ghf(1:end));                        
        gs = min([gs, 300]);
        ghx=hx(id(1:gs),:);  ghf=ghf(1:gs);           
        ghxd = real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
        spr = max(max(ghxd))/(Dim*gs)^(1/Dim);
        net = newrbe(ghx',ghf,spr);                       
        GlobalModelFUN = @(x) sim(net,x');  % build surrogate model                     
        % obtain candidate
        fitnessModel = GlobalModelFUN(U);                       
        [~,sidx] = sort(fitnessModel);                     
        candidate_position = U(sidx(1),:); 
        
    elseif  Action == 2     
        % Action 2: Surrogate sampling
        Action2 = Action2 + 1;
        [~, id] = sort(hf);
        lhx = hx(id(1:ls),:); lhf = hf(id(1:ls));
        lhxd = real(sqrt(lhx.^2*ones(size(lhx'))+ones(size(lhx))*(lhx').^2-2*lhx*(lhx')));
        spr = max(max(lhxd))/(Dim*ls)^(1/Dim);
        a1 = 100;
        net = newrbe(lhx',lhf,spr);       
        LocalModelFUN = @(x) sim(net,x');  % build surrogate model
        % obtain candidate
        Max_NFE = a1*Dim+1000; minerror = 1e-20;
        [candidate_position,~] = JADE(Dim, Max_NFE, LocalModelFUN, minerror, lhx); % find a optimum of surrogate model by optimizer
        
    elseif  Action == 3      
        % Action 3: Full-crossover 
        Action3 = Action3 + 1;
        [lhf, id] = sort(hf);               
        lhx = hx(id(1:ls2),:); lhf = hf(id(1:ls2));  
        lhxd = real(sqrt(lhx.^2*ones(size(lhx'))+ones(size(lhx))*(lhx').^2-2*lhx*(lhx')));
        spr = max(max(lhxd))/(Dim*ls2)^(1/Dim); 
        net = newrbe(lhx',lhf,spr);        
        LocalModelFUN = @(x) sim(net,x');  % build surrogate model
        % obtain candidate
        [candidate_position] = full_crossover(LocalModelFUN,lhx); 
        
    elseif  Action == 4      
        % Action 4: STR 
        Action4 = Action4 + 1;
        [lhf, id] = sort(hf);  
        lhx = hx(id(1:ls2),:); lhf = hf(id(1:ls2));  
        [newdata_x, newdata_f] = STR(FUN,lhx,lhf);
        % Update database and display 
        num_c = size(newdata_f,1);
        for a = 1:num_c
            NFE = NFE + 1; 
            candidate_position = newdata_x(a,:);
            candidate_fit = newdata_f(a);
            if show 
                disp(['  Fitness = ' num2str(candidate_fit) ' Solution = ' num2str(candidate_position) ]);
            end
            hx = [hx; candidate_position];  hf = [hf, candidate_fit];  
            [hf,idx] = sort(hf);                                        
            hx = hx(idx,:); 
            if  candidate_fit <= hf(1) 
                action_success = 1; R = 1;
                disp(['  Best fitness(Action ' num2str(Action)  ') = ' num2str(candidate_fit) ' NFE=' num2str(NFE)]);
            end
            Archive_FEs(NFE,:) = [NFE,candidate_fit];
            Archive_convergence(1,NFE) = min(Archive_FEs(1:NFE,2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update agent
        State_Next = 2*Action+action_success-1; 
        temp = max(Q_Agent(State_Next,:));
        Q_Agent(State,Action) = (1-alp)*Q_Agent(State, Action)+alp*(R+gamma*temp); 
        State = State_Next;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update agent
        continue;
        
    end
    
    % Update database and display 
    [~,ih,~] = intersect(hx,candidate_position,'rows');  
    if isempty(ih)~=1    % judge Repeat Sample 
        disp(['Sample repeat and delete it']);
        continue;
    end
    candidate_fit = FUN(candidate_position);
    NFE = NFE + 1; 
    if show 
        disp(['  Fitness = ' num2str(candidate_fit) ' Solution = ' num2str(candidate_position) ]);
    end
    hx = [hx; candidate_position];  hf = [hf, candidate_fit];   
    [hf,sidx] = sort(hf);                                     
    hx = hx(sidx,:); 
    if  candidate_fit <= hf(1) 
        action_success = 1; R = 1;
        disp(['  Best fitness(Action ' num2str(Action)  ') = ' num2str(candidate_fit) ' NFE=' num2str(NFE)]);
    end
    Archive_FEs(NFE,:) = [NFE,candidate_fit];
    Archive_convergence(1,NFE) = min(Archive_FEs(1:NFE,2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update agent
    State_Next = 2*Action+action_success-1; 
    temp = max(Q_Agent(State_Next,:));
    Q_Agent(State,Action) = (1-alp)*Q_Agent(State, Action)+alp*(R+gamma*temp); 
    State = State_Next;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update agent
end