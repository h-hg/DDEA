function EF_SAEA(Global)
% <algorithm> <A>
% A Ensemble Surrogate-based Framework for Expensive Multiobjective Optimization Problems
% dr --- 0.5 --- Auxiliary model

%------------------------------- Reference --------------------------------

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lay.Wu

%% Parameter setting
THETA0 = 5.*ones(Global.M,Global.D);
D = Global.D;
dr = floor(0.5 * D);
infill = zeros(5,D);
% np = 11*D-1; nv = D;
% THETA0 = (np^(-1/nv))*ones(Global.M, nv);
THETA.Theta1 = THETA0;
THETA.Theta2 = THETA0;
THETA.Theta3 = THETA0;
% THETA.Theta1 = THETA0(:,1:dr);
% THETA.Theta2 = THETA0(:,1:dr);

%% Generate the weight vectors and random population
[V0,Global.N] = UniformPoint(Global.N,Global.M);
N            = 11*Global.D-1;
PopDec       = lhsamp(N,Global.D);
Population   = INDIVIDUAL(repmat(Global.upper-Global.lower,N,1).*PopDec+repmat(Global.lower,N,1));
V = V0;
GDec = [];  Flag = 1;
FPopulation = Population;
%% Optimization
while Global.NotTermination(FPopulation)
    infill = zeros(5,D);
    PDec   = Population.decs;  PObj = Population.objs;
    [~,distinct1] = unique(round(PDec,6),'rows');
    PDec     = PDec(distinct1,:);  PObj = PObj(distinct1,:);
    PopDec = PDec;  PopObj = PObj;
    
    [Ensemble,r,THETA] = Ensemble_Construction(PopDec,PopObj,THETA,dr,Global);
    
%     [Dec,Obj] = Surrogate_assisted_Search1(PopDec,PopObj,Ensemble,V,r,Global);
%     [Dec,Obj] = Sur_NSGA2(PopDec,PopObj,Ensemble,r,Global);
%     [Dec,Obj] = Sur_MOEAD(PopDec,PopObj,Ensemble,r,Global,V);
    [Dec,Obj] = Sur_RVEA(PopDec,PopObj,Ensemble,r,Global,V);
%     [Dec,Obj] = Sur_NSGA3(PopDec,PopObj,Ensemble,r,Global);
    
    [infill,obj] = Infill_Solutions_Selection(Dec,Obj,PopDec,PopObj,V);
    
    infill1 = unique(infill,'rows');
    lia = ismember(infill1,Population.decs,'rows');
    infill = infill1(~lia,:);
    Isize = size(infill,1);
    if Isize > 0
        IPopulation = INDIVIDUAL(infill);
        InfObj = IPopulation.objs;
        Population = [Population,IPopulation];
    end
    
    [F1,~] = NDSort(Population.objs,1);
    FPopulation = Population(F1==1);
end
end