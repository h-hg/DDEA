function [ time,POP ] =RFCMOCO(problem_name)
% Usage: [ time,POP ] =RFCMOCO(problem_name)
% -----------------------------------------------------------------
% Input:
% problem_name  - File name of MOKP
%
% Output: 
% time          - Execution time
% POP           - Non-dominated solutions
%
    %%%%    Authors:    Handing Wang, Yaochu Jin
    %%%%    University of Surrey, UK
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage
    %%%%    DATE:       OCT 2018
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:
%Handing Wang, Yaochu Jin, A Random Forest-Assisted Evolutionary Algorithm %for Data-Driven Constrained Multiobjective Combinatorial Optimization of Trauma Systems IEEE Transactions on Cybernetics, accept, 2018.
%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------

rand('state',sum(100*clock));
load(problem_name);
c=size(p,2);%dimension of decision variable
ntree=100;% number of trees for random forest

pm=0.3;
n=100;
gmax=100;
NT=1000;
tic;
Bd  = ones(1,c); 
t=0.05;%significant level of LR
%Prepare the train data for RF
TrainData=Initialization( NT,c,Bd );
[ obj,Constraints ] =ComputeObjectivesConstraintsMOKP(TrainData,problem_name);
TrainData=[TrainData,obj,Constraints];
m=size(obj,2);
nc=size(Constraints,2);

[ Trees,S,NCtree] = SurrogateObjConRF(TrainData,c,m,nc,ntree);

I=randperm(size(TrainData,1));
[ Constraints] =PredictConstraintsMOKP(TrainData(I(1:50),1:c),Trees,S,NCtree,c,m,nc);
[ Beta,Derta ] =LR( TrainData(I(1:50),:),Constraints,c,t );

TPOP= Nondominated_C( TrainData,c,m,nc);
%Initialize the population 
POP = Initialization( n-size(TPOP,1),c,Bd );
POP=[POP;TPOP(:,1:c)];

g=1;
NG=5;
while g<=gmax&NT<=1500
    [ obj,Constraints] =PredictObjectivesConstraintsMOKP(POP(:,1:c),Trees,S,NCtree,c,m,nc);
    POP=[POP(:,1:c),obj,Constraints];
    %variation
    NPOP1= Cross( POP,TPOP,c,Bd );
    [obj,Constraints] =PredictObjectivesConstraintsMOKP(NPOP1,Trees,S,NCtree,c,m,nc);
    NPOP1=[NPOP1,obj,Constraints];
    NPOP2=MutationPopulation( TPOP,c,pm,Bd,n );
    [obj,Constraints] =PredictObjectivesConstraintsMOKP(NPOP2,Trees,S,NCtree,c,m,nc);
    NPOP2=[NPOP2,obj,Constraints];
    NPOP=[NPOP1;NPOP2];

    [ POP ] = POPcombination( POP,NPOP,c );
    
    %constraint correction
    ConLR=POP(:,c+m+1:c+m+nc)-ones(size(POP,1),1)*Derta;
    ConLR(find(ConLR<0))=0;
    POP(:,c+m+1:c+m+nc)=ConLR;
    %Stochastic Rank
    [ POP ] = FNDS( POP,c,m+nc);
    [ POP2 ] = FNDS_C( POP(:,1:c+m+nc),c,m,nc);
    [RS] = StochasticRanking( POP(:,1+c+m+nc),POP2(:,1+c+m+nc));
    NPOP=POP(RS,1:c+m+nc);
    
    %Model management
    [obj,s ] = PredictObjectives(TPOP(:,1:c),Trees,S,NCtree,c,m,nc);
    E = RMSE( obj,TPOP(:,c+1:c+m) );
    SPOP=NPOP(:,:);
    SPOP(:,c+1:c+m)=SPOP(:,c+1:c+m)-ones(size(SPOP,1),1)*E;
    if ~isempty(find(Derta<0))
        Derta=max(SPOP(:,c+m+1:c+m+nc),[],1)*0.1;
    end
    [ IR,ID] = Promissing_C( SPOP,TPOP,c,m,nc);
    SPOP1=SPOP(IR,:);
    SPOP([IR;ID],:)=[];
    [ TrainData] = POPcombination( TrainData,SPOP1,c );
    while size(TrainData,1)-NT<=NG & ~isempty(SPOP)
        ConLR=SPOP(:,c+m+1:c+m+nc)-ones(size(SPOP,1),1)*Derta*0.1;
        ConLR(find(ConLR<0))=0;
        SPOP(:,c+m+1:c+m+nc)=ConLR;
        [ IR,ID] = Promissing_C( SPOP,TPOP,c,m,nc);
        SPOP1=SPOP(IR,:);
        SPOP([IR;ID],:)=[];
        [ TrainData] = POPcombination( TrainData,SPOP1,c );
    end
    [TrainData(NT+1:end,c+1:c+m),TrainData(NT+1:end,c+m+1:c+m+nc)] =ComputeObjectivesConstraintsMOKP(TrainData(NT+1:end,1:c),problem_name);
    
    NT=size(TrainData,1);
    TPOP= Nondominated_C( TrainData,c,m,nc);
    %re-train RF
   [ Trees,S,NCtree] = SurrogateObjConRF(TrainData,c,m,nc,ntree);
    I=randperm(size(TrainData,1));
    [ Constraints] =PredictConstraintsMOKP(TrainData(I(1:50),1:c),Trees,S,NCtree,c,m,nc);
    [ Beta,Derta ] =LR( TrainData(I(1:50),:),Constraints,c,t );
    POP=NPOP(1:n,1:c+m+nc);
    
    g=g+1
end
[ POP(:,c+1:c+m),POP(:,c+m+1:c+m+nc) ] =ComputeObjectivesConstraintsMOKP(POP(:,1:c),problem_name);
[ POP ] = POPcombination( TrainData,POP,c );
POP= Nondominated_C( POP,c,m,nc);
toc;
time=toc;
end

