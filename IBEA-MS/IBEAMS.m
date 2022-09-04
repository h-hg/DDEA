function IBEAMS(Global)
    %------------------------------- Reference --------------------------------
    % Z. Liu, H. Wang, and Y. Jin, Performance Indicator based 
    % Adaptive Model Selection for Offline Data-Driven Multi-Objective 
    % Evolutionary Optimization in IEEE Transactions on Cybernetics.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2022 HandingWangXD Group. Permission is granted to copy and
    % use this code for research, noncommercial purposes, provided this
    % copyright notice is retained and the origin of the code is cited. The
    % code is provided "as is" and without any warranties, express or implied.
    %---------------------------- Parameter setting ---------------------------
    
    % IniN = 11*D-1-----The number of initial offline data
    % N    = 100--------The size of population
    % kappa= 0.05-------The parameter of IBEA parameter 
    % epsilon = 10^-5---The defined minimal distance difference and the sum
    % of RMSE
    % Objective relaxation: M-1
    % This code is written by Zhening Liu.
    % Email: zheningliu2@gmail.com
    
    %% Parameter setting
    Global.D = 20; 
    Global.M = 3;
    Global.N = 100;                     %Population size
    kappa    = 0.05;                    %The parameter of IBEA
    Generations    = 100;
    CurGen         = 1;
    Global.problem = 'DTLZ1';
    Global.lower   = zeros(1,Global.D);
    Global.upper   = ones(1,Global.D);
    cd(fileparts(mfilename('fullpath')));
    addpath(genpath(cd));
    %% Construct the kriging and RBFN models 
    PopDec      = LHS_sam(Global);                       %The initial offline data is sampled by the LHS
    Population  = Fitness(PopDec,Global);                %Generate the initial offline data evaulated by the real objective functions
    KModel      = construct_kriging(Population, Global);
    Rnets       = construct_Rnets(Population, Global);  
    MSE         = zeros(length(Population),Global.M);
    %% Optimization
    while CurGen < Generations
        MatingPool         = TournamentSelection(2,Global.N,-CalFitness(Population,kappa));
        OffDec             = GA(Global, decs(Population(MatingPool)), {1,20,1/Global.D,20});    
        [Population, MSE]  = AmendKriCal(Population, KModel, Global, MSE);                                  %Re-evaluate the parent population by the Kriging models
        [Offspring,OffMSE] = kriging_cal(OffDec, KModel, Global);                                           %Evaluate the Offspring by the Kriging models
        KFlag              = JudgeModel([Population,Offspring],[MSE;OffMSE]);                              %Select the Models to lead the optimization  
        if KFlag
            [Population,MSE] = EnvironmentalSelection([Population,Offspring],Global.N,kappa,[MSE;OffMSE]);  %Environmental selection assisted by the Kriging models
        else
            Population = AmendRBFCal(Population, Rnets, Global, MSE);                                       %Re-evaluate the parent population by the RBFN models
            Offspring  = RBFN_cal(OffDec, Rnets, Global);                                                %Re-evaluate the offspring by the RBFN models
            OffMSE     = ones(size(OffDec,1),Global.M);                                                     
            [Population,MSE] = EnvironmentalSelection([Population,Offspring],Global.N,kappa,[MSE;OffMSE]);  %Environmental selection assisted by the RBFN models
        end
        CurGen = CurGen+1;
    end
    OutPopulation = Population;
end
