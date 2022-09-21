 function [bestP,bestFitness] = PSO(Dim, Max_NFEs,FUN,minerror,ghx)
   
% disp('PSO global search');
n = Dim; NP = 100; flag_er=0;
lu = [min(ghx); max(ghx)];
LowerBound = lu(1, :);
UpperBound = lu(2, :);
 
%% PSO Parameters
VarSize = Dim;
MaxIt = 1000;
VarMax = UpperBound;
VarMin = LowerBound;
NFE = 0;

% Constriction Coefficients
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));
w=chi;          % Inertia Weight
wdamp=1;        % Inertia Weight Damping Ratio
c1=chi*phi1;    % Personal Learning Coefficient
c2=chi*phi2;    % Global Learning Coefficient

% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;

%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

particle=repmat(empty_particle,NP,1);

GlobalBest.Cost=inf;

for i=1:NP
    % Initialize Position
    particle(i).Position=(VarMax-VarMin).*rand(1,Dim)+VarMin;%随机产生初始种群个体
    
    % Initialize Velocity
    particle(i).Velocity=zeros(1,VarSize);
    
    % Evaluation
    particle(i).Cost=FUN(particle(i).Position);
    NFE = NFE +1;
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        GlobalBest=particle(i).Best;
    end
end
BestCost=zeros(MaxIt,1);
nfe=zeros(MaxIt,1);

%% PSO Main Loop
fitnessBest = GlobalBest.Cost;
time_begin=tic;
for it=1:MaxIt
    fitnessBest_old = fitnessBest;
    
    for i=1:NP
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(1,VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(1,VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
    end
    for i=1:NP
        Position(i,:) =  particle(i).Position;
    end
    
    cost = FUN(Position); NFE = NFE +NP;
    
    
    for i=1:NP
        particle(i).Cost =  cost(i);
    end
        
    for i=1:NP
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
        recRMSE(NFE)=particle(i).Cost;
    end
    
    
    BestCost(it)=GlobalBest.Cost;
    w=w*wdamp;
        fitnessBest = BestCost(it);
        error=abs(fitnessBest_old-fitnessBest); 
        if error <= minerror   
            flag_er=flag_er+1;
        else
            flag_er=0;
        end
        if flag_er >=10
            break;
        end
%     disp(['Iteration ' num2str(it) ' Best Cost = ' num2str(BestCost(it)) ' NFE=' num2str(NFE)]);
    bestFitness=BestCost(it);
    endNFEs = NFE;
    
    if NFE>Max_NFEs
       break
    end     
end
time_cost=toc(time_begin);
bestP=GlobalBest.Position;
