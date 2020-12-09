function [ time,gbest,POP ] =PSO_DV(Data,bu,bd )
% Usage: [ time,gbest,POP ] =PSO_DV(Data,bu,bd )
% -----------------------------------------------------------------
% Important Note: This code needs intralled SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:
% Data          - Data with c Decision Variables and Exact Objective Value
% bd            - Upper Boundary of c Decision Variables
% bd            - Lower Boundary of c Decision Variables
%
% Output: 
% time          - Execution Time
% Record        - Logbook of Evaluated Solutions
% gbest         - Predicted Optimum over Generations with c Decision Variables
%
    %%%%    Authors:    Handing Wang, Yaochu Jin, John Doherty
    %%%%    University of Surrey, UK
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage
    %%%%    DATE:       May 2018
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:
%Handing Wang, Yaochu Jin, John Doherty, Committee-based Active Learning for Surrogate-Assisted Particle Swarm Optimization of Expensive Problems, IEEE Transactions on Cybernetics, vol.47, no.9, pp.2664-2677, 2017.
%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
c=size(Data,2)-1;
x=Data(:,1:c);
y=Data(:,c+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit surrogates

% kriging
srgtOPTKRG  = srgtsKRGSetOptions(x, y);
srgtSRGTKRG = srgtsKRGFit(srgtOPTKRG);

% polynomial response surface
srgtOPTPRS  = srgtsPRSSetOptions(x, y);
srgtSRGTPRS = srgtsPRSFit(srgtOPTPRS);


% radial basis function
srgtOPTRBF  = srgtsRBFSetOptions(x, y);
srgtSRGTRBF = srgtsRBFFit(srgtOPTRBF);



n=50;

POP = initialize_pop(n,c,bu,bd);

[f1 PredVar] = srgtsKRGPredictor(POP(:,1:c), srgtSRGTKRG);
[f2 PredVar] = srgtsPRSPredictor(POP(:,1:c), x, srgtSRGTPRS);
f3    = srgtsRBFEvaluate(POP(:,1:c), srgtSRGTRBF);
F=[f1,f2,f3];
obj=(min(F,[],2)-max(F,[],2));%./mean(F,2);%0-std(F,[],2)./obj;%./mean(F,2);

POP=[POP,obj];
lbest=POP;

[best,Ib]=min(POP(:,c+1));
gbest=POP(Ib,:);
v=(2*rand(n,c)-1).*(ones(n,1)*(bu-bd))*0.5;
g=0;
gmax=100;
B=gbest;
while g<gmax  
    [ POP,v ] = fly(POP,bu,bd,gbest,lbest,v,g/gmax);

    
    [f1 PredVar] = srgtsKRGPredictor(POP(:,1:c), srgtSRGTKRG);
    [f2 PredVar] = srgtsPRSPredictor(POP(:,1:c), x, srgtSRGTPRS);
    f3    = srgtsRBFEvaluate(POP(:,1:c), srgtSRGTRBF);
    F=[f1,f2,f3];
    POP(:,c+1)=(min(F,[],2)-max(F,[],2));%./mean(F,2);%0-std(F,[],2)./obj;%./mean(F,2);
    [ POP,gbest,lbest] = Update_best(POP,gbest,lbest);
    best=gbest(end);
    if best<B(end,end)
        B=[B;gbest];
    end
    g=g+1;
end


time=0;

end