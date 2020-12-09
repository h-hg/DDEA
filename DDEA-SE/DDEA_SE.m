function [time,P,gbest ] =DDEA_SE(c,L,bu,bd)
% Usage: [time,P,gbest ] =DDEA_SE(c,L,bu,bd)
%
% Input:
% L             - Offline Data with c Decision Variables and Exact Objective Value
% c             - Number of Decision Variables
% bu            - Upper Boundary of c Decision Variables
% bd            - Lower Boundary of c Decision Variables
%
% Output: 
% time          - Execution Time
% P             - Final Predicted Optimum with c Decision Variables
% gbest         - Predicted Optimum over Generations with c Decision Variables
%
    %%%%    Authors:    Handing Wang, Yaochu Jin, Chaoli Sun, John Doherty
    %%%%    University of Surrey, UK and Taiyuan University of Science and Technology, China.
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage
    %%%%    DATE:       May 2018
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:

%Handing Wang, Yaochu Jin, Chaoli Sun, John Doherty, Offline data-driven evolutionary optimization using selective surrogate ensembles, IEEE Transactions on Evolutionary Computation, Accepted.

%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
rand('state',sum(100*clock));
nc=c;%Number of neurons of RBF models
T=2000;%Number of RBF models
Q=100;%Ensemble size
%Build Model Pool
[ W,B,C,S] = RBF_EnsembleUN( L,c,nc,T);
I=[1:T];%Using all the RBF Models
%-------------------------------------------
gmax=100;
pc=1;%Crossover Probability
pm=1/c;%Mutation Probability
n=100;%Population Sixe
%Online Optimization-------------------------------------------
tic;
POP = initialize_pop(n,c,bu,bd);
%RBF Predictors 
Y= RBF_Ensemble_predictor( W(I,:),B(I),C(:,:,I),S(:,I),POP,c );
POP=[POP,Y];
g=1;
gbest=[];
while g<=gmax
    %Model Management   
    if g~=1
        I= SelectModels(W,B,C,S,P(:,1:c),c,Q);
        POP=POP(:,1:c);
        Y= RBF_Ensemble_predictor(W(I,:),B(I),C(:,:,I),S(:,I),POP,c );
        POP=[POP,Y];
    end
    %Variations    
    NPOP1=SBX( POP,bu,bd,pc,n );
    [ Y ] = RBF_Ensemble_predictor( W(I,:),B(I),C(:,:,I),S(:,I),NPOP1,c );
    NPOP1=[NPOP1,Y];
    NPOP2=mutation(POP,bu,bd,pm,n);
    [ Y ] = RBF_Ensemble_predictor( W(I,:),B(I),C(:,:,I),S(:,I),NPOP2,c );
    NPOP2=[NPOP2,Y];
    POP=[POP;NPOP1;NPOP2];
    %Model Combination
    YAVE=mean(POP(:,c+1:end),2);
    [A,Is]=sort(YAVE);
    POP=[POP(Is(1:n),1:c)];
    g=g+1;
    P= POP(1,1:c);
    gbest=[gbest;P];

end

toc;
time=toc;

end

