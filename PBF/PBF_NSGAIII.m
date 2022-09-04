% Pareto-based Bi-indicator Infill sampling criterion based NSGAIII (PB-NSGAIII)
%------------------------------- Reference --------------------------------
% Z. Song, H. Wang, H. Xu, A Framework for Expensive Many-Objective
% Optimization with Pareto-based Bi-indicator Infill Sampling Criterion
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 HandingWangXD Group. Permission is granted to copy and
% use this code for research, noncommercial purposes, provided this
% copyright notice is retained and the origin of the code is cited. The
% code is provided "as is" and without any warranties, express or implied.

% This code is written by Zhenshou Song.
% Email: zssong@stu.xidian.edu.cn

clc;clear;warning off;

maxFes = 300;
popsize   = 100;
problem.problem = 'DTLZ1';
problem.M       = 3;
problem.D       = 10;
problem.lower   = zeros(1,problem.D);
problem.upper   = ones(1,problem.D);


%% PB-NSGAIII code
wmax = 15;

N = popsize;
M = problem.M;
P = lhsamp(N,problem.D);
Population = fitness (repmat(problem.lower, N, 1)+(repmat(problem.upper - problem.lower, N, 1)).*P, problem);

[W,~] = UniformPoint(N,problem.M);
Zmin  = min(Population.obj,[],1);
Training_data = Population;

THETA = 5.*ones(problem.M,problem.D);
Model = cell(1,problem.M);
k = 0;
Fes = N;
tic;
clc;fprintf('%s on %d-objective %d-variable (%6.2f%%), %.2fs passed...\n',problem.problem,problem.M,problem.D,Fes/maxFes*100,toc);

while Fes < maxFes
    
    Dec = Population.dec;
    Obj = Population.obj;
    train_X = Training_data.dec;
    train_Y = Training_data.obj;
    [~,distinct] = unique(roundn(train_X,-6),'rows');  
    train_X    = train_X(distinct,:);
    train_Y   = train_Y(distinct,:);
    for i = 1:problem.M % train surrogates
        dmodel     = dacefit(train_X,train_Y(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,problem.D),100.*ones(1,problem.D));
        Model{i}   = dmodel;
        THETA(i,:) = dmodel.theta;
    end
    
    w = 1;
    while w < wmax  % generate candidate population
        w  = w+1;
        MatingPool = TournamentSelection(2,popsize,zeros(popsize,1));
        OffspringDec  = GA(Dec(MatingPool,:),{1,20,1,20},problem);
        OffspringObj = zeros(size(OffspringDec,1),problem.M);
        for i = 1:size(OffspringDec,1)
            for j = 1:problem.M
                [OffspringObj(i,j),~,~] = predictor(OffspringDec(i,:),Model{j});
            end
        end
        
        Zmin  = min([Zmin; OffspringObj],[],1);
            
       %% Non-dominated sorting
        all_Obj = [Obj;OffspringObj];
        all_Dec = [Dec;OffspringDec];
        [FrontNo,MaxFNo] = NDSort([Obj;OffspringObj],popsize);
        Next = FrontNo < MaxFNo;    
        %% Select the solutions in the last front
        Last   = find(FrontNo==MaxFNo);  
        Choose = LastSelection(all_Obj(Next,:),all_Obj(Last,:),popsize-sum(Next),W,Zmin); 
        Next(Last(Choose)) = true;
        Dec = all_Dec(Next,:);
        Obj = all_Obj(Next,:);
    end   
    %%% PBISC
    DAdec = Dec;
    DA = Population;
    DA_Nor = (DA.obj - repmat(min([Obj;DA.obj],[],1),length(DA),1))...  % normalization
        ./repmat(max([Obj;DA.obj],[],1) - min([Obj;DA.obj],[],1),length(DA),1);
    DA_Nor_pre = (Obj - repmat(min([Obj;DA.obj],[],1),size(Obj,1),1))...
        ./repmat(max([Obj;DA.obj],[],1) - min([Obj;DA.obj],[],1),size(Obj,1),1);
    Zmin1 = min([DA_Nor;DA_Nor_pre],[],1);
    
    dist = zeros(size(DA_Nor_pre,1),size(DA_Nor,1));
    for i = 1:size(DA_Nor_pre,1)     % calculate the distance between candidate solutions and parents
        for j = 1:size(DA_Nor,1)
            dist(i,j) = norm(DA_Nor_pre(i,:)-DA_Nor(j,:),2);
        end
    end
    F1 = min(dist,[],2); % DI
    dist_D = pdist2(DA_Nor_pre,repmat(Zmin1,size(DA_Nor_pre,1),1));  % calculate the distance between candidate solutions and ideal point
    F2 = dist_D(:,1);  % % convergence indicator
    
    newObj= [(-F1),F2]; 
    ND1 = NDSort(newObj,1);
    PnewDec = DAdec((ND1==1),:);
    PnewDec = unique(PnewDec,'rows');
    
    Offspring = fitness(PnewDec, problem);
    Training_data.dec = [Training_data.dec;Offspring.dec];
    Training_data.obj = [Training_data.obj;Offspring.obj];
    ALL.dec       = [Population.dec;Offspring .dec];
    ALL.obj       = [Population.obj;Offspring .obj];
    
    Population = EnvironmentalSelection(ALL,popsize,W,Zmin);
    Fes = Fes + size(PnewDec,1);
    Zmin  = min(Population.obj,[],1);
    clc; fprintf('%s on %d-objective %d-variable (%6.2f%%), %.2fs passed...\n',problem.problem,problem.M,problem.D,Fes/maxFes*100,toc);
end

output = Training_data; % final output population;
PF = UniformPoint(10000,problem.M)/2;
IGD_value = IGD(output.obj,PF);


function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);

    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end
