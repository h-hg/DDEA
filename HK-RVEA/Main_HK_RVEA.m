%% this is the main file for the HK-RVEA algorithm, 
% Chugh, T., Allmendinger, R., Ojalehto, V., Miettinen, K., 
% Surrogate-assisted Evolutionary Biobjective Optimization for Objectives 
% with Non-uniform Latencies, in "Proceedings of the Genetic and 
% Evolutionary Computation Conference (GECCO-2018)", 
% Ed. by H. Aguirre, The Association of Computing Machinery (ACM), 
% 609-616, 2018

%% 
% please see the licence file for more information
%% contact: Tinkle Chugh, Kaisa Miettinen, Vesa Ojalehto and Richard Allmendinger if you have any questions
% Emails: t.chugh@exeter.ac.uk, kaisa.miettinen@jyu.fi, richard.allmendinger@manchester.ac.uk

%% 
clear; clc; close all;
Problem = 'DTLZ2'; 
% the implenetation is for two objectives
id_ex = 2; % objective number for the most expensive obj function
id_nex = 1;
addpath(genpath('support_files')); 
%% Input: 
M = 2; % number of objectives
no_var = 10; Bounds = [ones(1,no_var);zeros(1,no_var)]; % number of variables and their bounds
latency = 2; % Latency value
Max_FE_nex = 300; % 
Max_FE_ex = round(Max_FE_nex/latency);% maximum no of expensive evaluations for the most expensive objective function
%% Step-1: Generate the initial data

P = generate_initial_data(Bounds);
P = unique(P,'rows');

FE_ex = 0; FE_nex = 0; A = []; A_nex = []; A_ex = [];
itr_count = 1; empty_ref = 0;

while (FE_ex < Max_FE_ex || FE_nex < Max_FE_nex)
    
   
    %% step 2: evaluate P on the most expensive objective function - 
    F_exp = evaluate_most_expensive_obj(P,Problem,id_ex);
    A_ex = [A_ex;[P,F_exp]];
    
    %% step 3: 
 
    if latency>1
        if itr_count ==1 
            [X_nex,F_nex] = optimize_least_expensive(P,Bounds,latency,Problem,id_nex);   
        else
            [X_nex,F_nex] = genetic_operation(P,Bounds,latency,Problem,id_nex);
        end
    else
        
        F_nex = evaluate_least_expensive_obj(P,Problem,id_nex);
        X_nex = P;
    end

    A_nex = [A_nex;[X_nex,F_nex]];
    %% step 4:
    FE_ex = FE_ex + size(P,1);
    FE_nex = FE_nex + size(P,1)*latency
    A = [A;select_solutions_for_archive(P,F_exp,F_nex,id_ex,id_nex)];
    
    %% Step 5: Build surrogate for both objective functions
    
    model_ex = build_model(A_ex,no_var);
    model_nex = build_model(A_nex,no_var);
  
    %% Step 6: Run some surrogate-assisted algorithm e.g. K-RVEA to find the samples to be evaluated

    P = run_K_RVEA(model_ex,model_nex,Bounds,A,id_ex,id_nex,empty_ref);

    itr_count = itr_count + 1;
 
end

%% Output:
non = P_sort(A(:,no_var+1:end),'first')==1;
PF_A = A(non,:);
figure;
scatter(A(:,no_var+1),A(:,no_var+2));
hold on;
scatter(PF_A(:,no_var+1),PF_A(:,no_var+2));
hold off;
legend('All solutions','Nondominated solutions');

%% Support-files
function X = generate_initial_data(Bounds)
    no_var = size(Bounds,2);
    sample_size = 10*no_var;
    Xn = lhsdesign(sample_size,no_var);
    ub = Bounds(1,:); lb = Bounds(2,:);
    X = bsxfun(@plus,lb,bsxfun(@times,Xn,(ub-lb)));
end

function A = select_solutions_for_archive(P,F_exp,F_nex,id_ex,id_nex)
    S = size(P,1);  
    F = zeros(S,2);
    F(:,id_ex) = F_exp;
    F(:,id_nex) = F_nex(1:S,:);
  
    A = [P,F];
end

function model = build_model(A,no_var)
    
X_train = A(:,1:no_var);
Y_train = A(:,no_var+1:end);

    if size(A,1)>500
        X_temp = [];
        Y_temp = [];
        idx = kmeans(X_train,500);
        for n = 1:500
            t = find(idx==n); pos = randi(length(t));
            ind = t(pos);
            X_temp = [X_temp; X_train(ind,:)];    
            Y_temp = [Y_temp; Y_train(ind,:)];   
        end
    X_train = X_temp;
    Y_train = Y_temp;
    end
    model = fitrgp(X_train,Y_train,'KernelFunction','ardsquaredexponential');
end

function [pop, Fitness]  = genetic_operation(Population,Bounds,latency,Problem,id_nex)
        N = size(Population,1)*latency - size(Population,1);
        MatingPool= F_mating(Population,N);
        Coding = 'Real';
        Offspring = P_generator(MatingPool,Bounds,Coding,N);
        
        pop = [Population;Offspring];
        
        F = P_objective('value',Problem,2,pop);
        Fitness = F(:,id_nex);
end


