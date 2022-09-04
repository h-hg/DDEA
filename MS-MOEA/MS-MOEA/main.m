% An Adaptive Model Switch-based Surrogate-Assisted Evolutionary Algorithm 
% for Noisy Expensive Multi-Objective Optimization
%------------------------------- Reference --------------------------------
% N. Zheng, H. Wang, and B. Yuan, An Adaptive Model Switch-based 
% Surrogate-Assisted Evolutionary Algorithm for Noisy Expensive 
% Multi-Objective Optimization, Complex & Intelligent Systems, accepted, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 HandingWangXD Group. Permission is granted to copy and
% use this code for research, noncommercial purposes, provided this
% copyright notice is retained and the origin of the code is cited. The
% code is provided "as is" and without any warranties, express or implied.

% This code is written by Nan Zheng.
% Email: nanszheng@163.com

clc;
% clear all;
close all;
addpath(genpath(pwd));
%%  Parameters Setting
wmax=20;            % Maximum number of model-assisted iterations
mu=5;               % Maximum number of samples
kappa=0.05;         % Scaling factor
jia=1;              % Noise Sign 1/Within°¢0/Without
threshold=0.03;     % Noise threshold
fun_test='DTLZ4M2'; % Test function
choose_Nflag=0;     % 0/RBFN£¨1/Kriging
Offspring_lastc=[]; % Last convergence candidates
Offspring_lastd=[]; % Last diversity candidates
%% Select test function
[num_vari,num_obj,Size,max_evaluation,decision_space,PF]=benchmark_chosen(fun_test);
% Number of variables°¢Number of objectives°¢Population size°¢Maximum number of evaluations°¢Decision space°¢Pareto front%
THETA = 5.*ones(num_obj,num_vari);
%% Create initial population
[initial_x,initial_y]=Pop_initial(fun_test,Size,num_vari,num_obj,jia);
%% Individual archiving with noise true evaluation
T_Aori=[initial_x,initial_y];
for i=1:size(T_Aori,1)
    T_AIori(1,i)=INDIVIDUAL(T_Aori(i,1:num_vari),T_Aori(i,num_vari+1:end));
end

%% Estimating noise impact 
           test_solution_obj=zeros(5,num_obj);
           test_solution=Evolutionaryperiod(T_AIori,kappa,45);
           test_solution_dec=test_solution.dec;
           test_solution_obj(1,:)=test_solution.obj;           
           for i=1:4
               if(num_obj==2)
                  test_solution_obj(i+1,:) = feval(fun_test, test_solution_dec, num_obj)+mvnrnd ([0 0],[0.2 0 ;0 0.2 ],size(test_solution_dec,1));
               else
                  test_solution_obj(i+1,:) = feval(fun_test, test_solution_dec, num_obj)+mvnrnd ([0 0 0],[0.2 0 0;0 0.2 0;0 0 0.2],size(test_solution_dec,1));
               end               
           end
           test_solution_obj_ave=sum(test_solution_obj)/size(test_solution_obj,1);
           test_solution_obj_err=abs(test_solution_obj-repmat(test_solution_obj_ave,size(test_solution_obj,1),1));
           test_solution_max_err=max(test_solution_obj_err);
           SNR_max=test_solution_max_err./abs(test_solution_obj_ave);
           for j=1:size(SNR_max,2)
               if SNR_max(1,j)<=threshold
                  choose_Nflag=1;break;
               end
           end
%% Data pre-processing
         [y_deal,T_A,T_A_RBF,T_A_Clu]=data_deal(T_Aori,T_AIori,num_obj,num_vari,choose_Nflag);
            
         for i=1:Size
            Population(1,i)=INDIVIDUAL(initial_x(i,:),y_deal(i,:));
         end
         for i=1:size(T_A,1)
             T_AI(1,i)=INDIVIDUAL(T_A(i,1:num_vari),T_A(i,num_vari+1:end));
         end
         evaluation=size(T_Aori,1)+4;
         PopulationDec=Population.decs;
         PopulationObj=Population.objs;
%% Optimization
while evaluation<max_evaluation  
        %% Noise impact judgment    
         if choose_Nflag==1            
              choose_Nflag=0;
              test_solution=Evolutionaryperiod(T_AI,kappa,Size);
              SNR_max=test_solution_max_err./abs(test_solution.obj);
              for j=1:size(SNR_max,2)
               if SNR_max(1,j)<=threshold
                  choose_Nflag=1;break;
               end
              end
         end
         %% Construct surrogate model 0/RBFN£¨1/Kriging
          y_deal=[];
         if choose_Nflag==1
            for i = 1 : num_obj
                dmodel     = dacefit(T_A(:,1:num_vari),T_A(:,num_vari+i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,num_vari),100.*ones(1,num_vari));
                model_K{i}   = dmodel;
                THETA(i,:) = dmodel.theta;
            end            
         end 
         % surrogate model and smoothing model         
          for i = 1 : num_obj
            dmodel     = construct_rbfn(T_Aori(:,1:num_vari),T_Aori(:,num_vari+i),evaluation);
            Model{i}   = dmodel;
          end
                    
          if choose_Nflag==1
            for i = 1: size(PopulationDec,1)
                for j = 1 : num_obj
                    [PopulationObj(i,j),~,Popmse(i,j)] = predictor(PopulationDec(i,:),model_K{j});
                end
            end
          end
        if choose_Nflag==0 
          for j = 1 : num_obj
            PopulationObj(:,j) = cal_via_net(PopulationDec, Model{j});
          end
        end
         for i=1:Size
            Population(1,i)=INDIVIDUAL(PopulationDec(i,:),PopulationObj(i,:));
         end
      %% Model-assisted evolutionary search  
      w  = 1;
      while w <= wmax
      if w==1
      MatingPool = TournamentSelection(2,Size,-CalFitness(PopulationObj,kappa));
      else    
      MatingPool = TournamentSelection(2,size(Offspring,2),-CalFitness(OffspringObj,kappa));
      end
      if w==1 
       OffspringDec  = GA(Population(MatingPool),num_vari); 
      else
       OffspringDec  = GA(Offspring(MatingPool),num_vari);
      end
      [N,~]  = size(OffspringDec);
      OffspringObj = zeros(N,num_obj); 
      %% RBFN as surrogate model
       if choose_Nflag==0 
            for j = 1 : num_obj
                OffspringObj(:,j) = cal_via_net(OffspringDec, Model{j});
            end
            if w == 1
                Offspring_Old=Population;                
            else
                Offspring_Old=Offspring;
            end
            for i=1:size(OffspringDec,1)
                Offspring(1,i)=INDIVIDUAL(OffspringDec(i,:),OffspringObj(i,:));
            end
            Offspring_All=[Offspring_Old,Offspring];
            if w==wmax
               Offspring_candidate=Offspring_All;
            end
            % RBFN as surrogate model,Selection based on diversity
            Offspring1 = RModelEnvironmentalSelectionD(Offspring_All,Size-50);           
            for i=1:size(Offspring1,2)
                Offspring_All(find(Offspring_All==Offspring1(i)))=[];
            end
            % RBFN as surrogate model,Selection based on convergence
            Offspring2 = RModelEnvironmentalSelectionC(Offspring_All,Size-size(Offspring1,2),kappa);             
            Offspring=[Offspring1,Offspring2];
            OffspringObj=Offspring.objs;
            OffspringDec=Offspring.decs;    
       end      
      %% Kriging as surrogate model
      if choose_Nflag==1
            if w == 1
                Offspring_Old=Population;
                MSE_Old=Popmse;
            else
                Offspring_Old=Offspring;
                MSE_Old=Offmse;
            end   
           for i = 1: size(OffspringDec,1)
                for j = 1 : num_obj
                    [OffspringObj(i,j),~,Offmse(i,j)] = predictor(OffspringDec(i,:),model_K{j});
                end
           end
            for i=1:size(OffspringDec,1)
                Offspring(1,i)=INDIVIDUAL(OffspringDec(i,:),OffspringObj(i,:));
            end
            Offspring_All=[Offspring_Old,Offspring];
            MSE_All=[MSE_Old;Offmse];
            if w==wmax
               Offspring_candidate=Offspring_All;
               Mse_candidate=MSE_All;
            end
            % Kriging as surrogate model,Selection based on convergence
            [Offspring1,Offmse1]=KModelEnvironmentalSelectionC(Offspring_All,MSE_All,0.5*Size,kappa);
            for i=1:size(Offspring1,2)
                k=find(Offspring_All==Offspring1(1,i));
                Offspring_All(k)=[];
                MSE_All(k,:)=[];
            end
            % Kriging as surrogate model,Selection based on uncertainty
            [Offspring2,Offmse2]= Selection_uncertainty(Offspring_All,MSE_All,0.5*Size);

            Offspring=[Offspring1,Offspring2];
            Offmse=[Offmse1;Offmse2];

            OffspringObj=Offspring.objs;
            OffspringDec=Offspring.decs;                  
      end
            w = w + 1;

      end
    %% RBFN as surrogate model,Selecting candidates based on convergence and diversity
    if  choose_Nflag==0
      Offspring_Covergence= RModelEnvironmentalSelectionC(Offspring_candidate,Size,kappa);   
      Offspring_Diverity=RModelEnvironmentalSelectionD(Offspring_candidate,Size);
    end
    
    %% Kriging as surrogate model,Selecting candidates based on convergence and diversity
    if  choose_Nflag==1
      [Offspring_Covergence,MSE_Covergence]=KModelEnvironmentalSelectionC(Offspring_candidate,Mse_candidate,Size,kappa);
      [Offspring_Diverity,MSE_Diverity]=KModelEnvironmentalSelectionD(Offspring_candidate,Mse_candidate,Size);   
    end
    
      %% Select the solution to be evaluated
    if choose_Nflag==0
       [PopNew,Offspring_lastc,Offspring_lastd,~]=SelectionN(Offspring_Covergence,Offspring_Diverity,kappa,Offspring_lastc,Offspring_lastd,Model,choose_Nflag);
    end
    if choose_Nflag==1
       [PopNew,Offspring_lastc,Offspring_lastd,~]=Selection(Offspring_Covergence,Offspring_Diverity,MSE_Diverity,kappa,Offspring_lastc,Offspring_lastd,model_K,choose_Nflag);
    end
    %% Delete the selected duplicate solutions to be evaluated
      PopNew1=SelectionNrepeat(PopNew);
      PopNew=[];PopNew=PopNew1;      
      PopNew_x=[];
    %% Delete solutions to be evaluated that are too close to each other
    for i = 1:size(PopNew,1)
        dist2 = pdist2(real( PopNew(i,:)),real(T_AI.decs));
        if min(dist2) > 1e-5
            PopNew_x = [PopNew_x;PopNew(i,:)];
        end
    end
    %% Function Evaluation
    if size(PopNew_x,1)~=0        
        if jia==1
           if (num_obj==2)
            PopNew_y  = feval(fun_test, PopNew_x, num_obj)+mvnrnd ([0 0],[0.2 0 ;0 0.2 ],size(PopNew_x,1));
           else
            PopNew_y  = feval(fun_test, PopNew_x, num_obj)+mvnrnd ([0 0 0],[0.2 0 0;0 0.2 0;0 0 0.2],size(PopNew_x,1));    
           end
           
        else
        PopNew_y  = feval(fun_test, PopNew_x, num_obj);
        end
    

               
       New       = [PopNew_x,PopNew_y];
       T_Aori        = [T_Aori;New];
        for i=1:size(T_Aori,1)
            T_AIori(1,i)=INDIVIDUAL(T_Aori(i,1:num_vari),T_Aori(i,num_vari+1:end));
        end

         %% Noise reduction treatment
         y_RBFdeal=[];y_Cludeal=[];
        for j = 1 : num_obj
            y_RBFdeal(:,j) = cal_via_net(T_Aori(:,1:num_vari), Model{j}); 
        end
            y_RBFdeal=y_RBFdeal;
            T_A_RBF        = [T_Aori(:,1:num_vari),y_RBFdeal]; 
            y_RBFdeal=T_A_RBF(:,num_vari+1:end);
            y_Cludeal=T_AIori.objs;
            T_A_Clu=[T_Aori(:,1:num_vari),y_Cludeal];
            if choose_Nflag==1
                y_deal=y_Cludeal;
                T_A=T_A_Clu;
            else
                y_deal=y_RBFdeal;
                T_A=T_A_RBF;
            end                    
            evaluation=size(T_Aori,1)+4;

       for i=1:size(T_A,1)
            T_AI(1,i)=INDIVIDUAL(T_A(i,1:num_vari),T_A(i,num_vari+1:end));
       end

        %% Selecting iterative populations
        PopulationDec=Selection_Pop(PopNew_x,PopNew_y,T_AI,Size,kappa);
        
     end        
    
end
%%  ‰≥ˆ÷÷»∫
OutPut=T_AI;






    
    