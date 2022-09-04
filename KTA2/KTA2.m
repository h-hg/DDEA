% Kriging assisted Two_Arch2
%------------------------------- Reference --------------------------------
% Z. Song, H. Wang, C. He and Y. Jin, A Kriging-Assisted Two-Archive
% Evolutionary Algorithm for Expensive Many-Objective Optimization in IEEE
% Transcations on Evolutionary Computation.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 HandingWangXD Group. Permission is granted to copy and
% use this code for research, noncommercial purposes, provided this
% copyright notice is retained and the origin of the code is cited. The
% code is provided "as is" and without any warranties, express or implied.
%---------------------------- Parameter setting ---------------------------
% mu   ---  5 --- Number of re-evaluated solutions at each generation
% tau  ---  0.75*N --- The proportion of one type noninfluential points in
% training data
% wmax ---  10 ---  Number of generations before updating CA and DA
% phi  ---  0.1 --- Number of randomly selected individuals

% This code is written by Zhenshou Song.
% Email: zssong@stu.xidian.edu.cn

clc; clear; warning off;



maxFes    = 300;
popsize   = 100;
problem.problem = 'DTLZ1';
problem.M       = 3;
problem.D       = 10;
problem.lower   = zeros(1,problem.D);
problem.upper   = ones(1,problem.D);



tau  = 0.75;
phi  = 0.1;
wmax = 10;
mu   = 5;
p = 1/problem.M;
CAsize = popsize;
N = popsize;
P = lhsamp(N,problem.D);
Population = fitness (repmat(problem.lower, N, 1)+(repmat(problem.upper - problem.lower, N, 1)).*P, problem);

Training_data = Population;
CA = UpdateCA([],Population,CAsize);
DA = Population;
THETA_S = 5.*ones(problem.M,problem.D);
THETA_IS =  5.*ones(2,problem.M,problem.D);
Model_sensitive = cell(1,problem.M);
Model_insensitive = cell(2,problem.M);
Fes = popsize;
tic;
clc; fprintf('%s on %d-objective %d-variable (%6.2f%%), %.2fs passed...\n',problem.problem,problem.M,problem.D,Fes/maxFes*100,toc);
while Fes < maxFes
    for i = 1:problem.M
        dmodel     = dacefit(Training_data.dec,Training_data.obj(:,i),'regpoly0','corrgauss',THETA_S(i,:),1e-5.*ones(1,problem.D),100.*ones(1,problem.D));
        Model_sensitive{i}   = dmodel;
        THETA_S(i,:) = dmodel.theta;
    end
    Centers = zeros(problem.M,2);
    for i = 1:problem.M
        [~,N1] = sort(Training_data.obj(:,i));
        num = ceil(size(Training_data.dec,1).*tau);
        mean_index{1} = N1(1:num);
        mean_index{2} = N1(end-num:end);
        for j = 1:2
            Centers(i,j) = mean(Training_data.obj(mean_index{j},i));  % lambda and miu
        end
        for j = 1:2
            train_X = Training_data.dec(mean_index{j},:);
            train_Y = Training_data.obj(mean_index{j},i);
            dmodel  = dacefit(train_X,train_Y,'regpoly0','corrgauss',THETA_IS(j,i,:),1e-5.*ones(1,problem.D),100.*ones(1,problem.D));
            Model_insensitive{j,i} = dmodel;
            THETA_IS(j,i,:)        = dmodel.theta;
        end
    end
    
    CAobj = CA.obj; CAdec = CA.dec;
    DAobj = DA.obj; DAdec = DA.dec;
    w = 1;
    while w <= wmax
        [~,ParentCdec,~,ParentMdec] = MatingSelection(CAobj,CAdec,DAobj,DAdec,popsize);
        OffspringDec = [GA(ParentCdec,{1,20,0,0},problem);GA(ParentMdec,{0,0,1,20},problem)];
        PopDec = [DAdec;CAdec;OffspringDec];
        N2      = size(PopDec,1);
        PopObj = zeros(N2,problem.M);
        MSE    = zeros(N2,problem.M);
        for i = 1:N2
            for j = 1:problem.M
                [PopObj(i,j),~,~] = predictor(PopDec(i,:),Model_sensitive{j});
                if abs(PopObj(i,j)- Centers(j,1)) <= abs(PopObj(i,j)- Centers(j,2))
                    model = Model_insensitive{1,j};
                else
                    model = Model_insensitive{2,j};
                end
                [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),model);
            end
        end
        [CAobj,CAdec,~] = K_UpdateCA(PopObj,PopDec,MSE,CAsize);
        [DAobj,DAdec,DAvar] = K_UpdateDA(PopObj,PopDec,MSE,popsize,p);
        w = w + 1;
    end
    
    Offspring01 = Adaptive_sampling(CAobj,DAobj,CAdec,DAdec,DAvar,DA,mu,p,phi);
    
    [~,index] = unique(Offspring01 ,'rows');
    PopNew = Offspring01(index,:);
    Offspring02 = [];
    for i = 1:size(PopNew,1)
        dist2 = pdist2(real( PopNew(i,:)),real(Training_data.dec));
        if min(dist2) > 1e-5
            Offspring02 = [Offspring02;PopNew(i,:)];
        end
    end
    
    if ~isempty(Offspring02)
        Offspring = fitness(Offspring02,problem);
        Fes = Fes + size(Offspring02,1);
        Training_data.dec = [Training_data.dec;Offspring.dec];
        Training_data.obj = [Training_data.obj;Offspring.obj];
        CA = UpdateCA(CA,Offspring,CAsize);
        DA = UpdateDA(DA,Offspring,popsize,p);
        clc; fprintf('%s on %d-objective %d-variable (%6.2f%%), %.2fs passed...\n',problem.problem,problem.M,problem.D,Fes/maxFes*100,toc);
    end
    
end

output = Training_data; % final output population;












