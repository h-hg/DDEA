%% This is the MATLAB code for the K-RVEA algorithm published in the following article:
%  T. Chugh, Y. Jin, K. Miettinen, J. Hakanen, and K. Sindhya, “A
% surrogate-assisted reference vector guided evolutionary algorithm for
% computationally expensive many-objective optimization,” IEEE Transactions 
% on Evolutionary Computation, vol. 22, no. 1, pp. 129–142, 2018.

% More details about it can be found in the thesis: T. Chugh, “Handling 
% expensive multiobjective optimization problems with evolutionary 
% algorithms,” Ph.D. dissertation, Jyväskylä Studies in
% Computing 263, University of Jyväskylä, 2017

% Please read the licence file before using the code and cite the article
% and the thesis if you use the code. 


% Contact: t.chugh@exeter.ac.uk, tinkle.chugh@gmail.com,
% kaisa.miettinen@jyu.fi
%%
clear;
clc;

addpath(genpath('Public'));
addpath(genpath('dace')); % library for Kriging
clear dmodel;
warning off;
Problems = {'DTLZ2'}; % problems to be solved
obj = [2,3,4,5,6]; % number of objectives

maxFE = 250; % maximum number of expensive function evaluations
maxrun = 21; % number of independent runs

Int_Col = 0;
Empty_ref_old = 0;
%%
for p = 1:length(Problems)
    Problem = Problems{p};      
    for run = 1:maxrun
        for k = 1:length(obj)
            M = obj(k);
            Vs = ref_vectors(M,Problem); % Generation of reference vectors
            [~,Boundary,Coding] = P_objective('init',Problem,M,size(Vs,1));
            no_var = size(Boundary,2);
            W = Vs;
           
            %% Generating the data
            sample_size = 11*no_var - 1;
            Xn = lhsdesign(sample_size,no_var);
            ub = Boundary(1,:);
            lb = Boundary(2,:);
            Population = bsxfun(@plus,lb,bsxfun(@times,Xn,(ub-lb)));
            FunctionValue = P_objective('value',Problem,M,Population);
            %% Main loop
            FE = size(Population,1);
            while FE <=maxFE
                [~,distinct] = unique(Population,'rows');
                Population = Population(distinct,:);
                FunctionValue = FunctionValue(distinct,:);
                % Dont change the parameters here
                theta1 = 0.5*ones(1,no_var); lob = 0.00001*ones(1,no_var); upb = 200*ones(1,no_var);
                for t=1:M
                    dmodel(t) = dacefit(Population,FunctionValue(:,t), @regpoly0, @corrgauss, theta1, lob, upb);
                end
                [pop,current_pop,Empty_ref_old] = evolve_K_RVEA(dmodel,Population,FunctionValue,Boundary,Int_Col,Vs,Empty_ref_old);
                Lia = ismember(pop,Population,'rows');
                r_unique = find(Lia(:,1)==0);
                pop = pop(r_unique,:);
                tt = 1;
                while (isempty(pop))
                    rt = randperm(size(current_pop,1));
                    pop = current_pop (rt(1:1),:);
                    Lia = ismember(pop,Population,'rows');
                    r_unique = find(Lia(:,1)==0);
                    pop = pop(r_unique,:);
%                     size (off)
                    tt = tt+1;
                    if tt>100
                        break;
                    end
                end    
                if tt<=100
                    fitness = P_objective('value',Problem,M,pop);
                    
                    FE = FE + size(pop,1)
                    FunctionValue = [FunctionValue;fitness];
                    Population = [Population;pop];
                end
            end
%             plot(t_record)
            solutions(run,k).(Problem) = [Population,FunctionValue]; 
        
        clear Population;
        clear FunctionValue;
        end
    end
    save ('KRVEA_solutions_so_far.mat','solutions')
end
save ('KRVEA_solutions_final.mat','solutions')
% toc;
