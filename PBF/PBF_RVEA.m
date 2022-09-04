% Pareto-based Bi-indicator Infill sampling criterion based RVEA (PB-RVEA)
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


%% PB-RVEA code
alpha = 2;
wmax = 15;

[V0,~] = UniformPoint(popsize,problem.M);
V1 = V0;
V = V1;
N = popsize;
M = problem.M;
P = lhsamp(N,problem.D);
Population = fitness (repmat(problem.lower, N, 1)+(repmat(problem.upper - problem.lower, N, 1)).*P, problem);


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
    while w <= wmax   % generate candidate population
        theta = (w/wmax)^alpha; 
        OffspringDec = GA(Dec,{1,20,1,20},problem);
        OffspringObj = zeros(size(OffspringDec,1),problem.M);
        for i = 1:size(OffspringDec,1)
            for j = 1:problem.M
                [OffspringObj(i,j),~,~] = predictor(OffspringDec(i,:),Model{j});
            end
        end
        
        PopObj = [Obj;OffspringObj];
        PopDec = [Dec;OffspringDec];
        PopObj1 = PopObj;
        
        [N,M]  = size(PopObj);
        NV     = size(V,1);
        PopObj = PopObj - repmat(min(PopObj,[],1),N,1);
        cosine = 1 - pdist2(V,V,'cosine');
        cosine(logical(eye(length(cosine)))) = 0;
        gamma  = min(acos(cosine),[],2);  
        %% Associate each solution to a reference vector
        Angle = acos(1-pdist2(PopObj,V,'cosine'));
        [~,associate] = min(Angle,[],2);
        Next = zeros(1,NV);
        
        for i = unique(associate)'
            current1 = find(associate==i);
            APD = (1+M*theta*Angle(current1,i)/gamma(i)).*sqrt(sum(PopObj(current1,:).^2,2));
            [~,best] = min(APD);
            Next(i)  = current1(best);     
        end
        index = Next(Next~=0);
        Dec = PopDec(index,:);
        Obj = PopObj1(index,:);
        if ~mod(w,ceil(wmax*0.1))
            V(1:size(V0,1),:) = V0.*repmat(max(Obj,[],1)-min(Obj,[],1),size(V0,1),1);
        end
        w = w +1;
    end   
    %%% PBISC
    DAdec = Dec;
    DA = Population;
    DA_Nor = (DA.obj - repmat(min([Obj;DA.obj],[],1),length(DA),1))...  % normalization
        ./repmat(max([Obj;DA.obj],[],1) - min([Obj;DA.obj],[],1),length(DA),1);
    DA_Nor_pre = (Obj - repmat(min([Obj;DA.obj],[],1),size(Obj,1),1))...
        ./repmat(max([Obj;DA.obj],[],1) - min([Obj;DA.obj],[],1),size(Obj,1),1);
    Zmin = min([DA_Nor;DA_Nor_pre],[],1);
    
    dist = zeros(size(DA_Nor_pre,1),size(DA_Nor,1));
    for i = 1:size(DA_Nor_pre,1)    % calculate the distance between candidate solutions and parents
        for j = 1:size(DA_Nor,1)
            dist(i,j) = norm(DA_Nor_pre(i,:)-DA_Nor(j,:),2);
        end
    end
    
    F1 = min(dist,[],2); % diversity indicator
    
    dist_D = pdist2(DA_Nor_pre,repmat(Zmin,size(DA_Nor_pre,1),1)); % calculate the distance between candidate solutions and ideal point
    
    F2 = dist_D(:,1);  % convergence indicator
    
    newObj= [(-F1),F2]; 
    ND1 = NDSort(newObj,1);
    PnewDec = DAdec((ND1==1),:);  % find solutions in the first front
    PnewDec = unique(PnewDec,'rows');
    
    Offspring = fitness(PnewDec, problem); % real function evaluation
    Training_data.dec = [Training_data.dec;Offspring.dec];
    Training_data.obj = [Training_data.obj;Offspring.obj];
    mu = size(PnewDec,1);
    NI = size(Population.dec,1);
    Population = UpdataArchive(Population,Offspring,V1,mu,NI); 
    V1(1:size(V0,1),:) = ReferenceVectorAdaptation(Population.obj,V0);
    k = k+1;
    V = V1; 
    Fes = Fes + size(PnewDec,1);
    clc; fprintf('%s on %d-objective %d-variable (%6.2f%%), %.2fs passed...\n',problem.problem,problem.M,problem.D,Fes/maxFes*100,toc);
end

output = Training_data; % final output population;
PF = UniformPoint(10000,problem.M)/2;
IGD_value = IGD(output.obj,PF);
