
IP_train = pastpos;
OT_train = pastposfit;
INput=[ OT_train IP_train];
PREIN=INput;
[PRETrainL,PRETrainD]=size(PREIN); 


st=0.8;
STId=rand(PRETrainL,1);%STId is the id of slected train data
I1=find(STId<=st);
I0=find(STId>st);

P_train =PREIN(I1,(2:D+1));
T_train =PREIN(I1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_test = p;
sigma0 = std(T_train);




gprMdl = fitrgp(P_train,T_train,'KernelFunction','matern32','BasisFunction','pureQuadratic','FitMethod','sr', 'Standardize',1, 'ComputationMethod','v');


[t_sim,mse]=predict(gprMdl,T_test);

posfit=t_sim;
posmse=mse;
 PopObj=[t_sim mse];

fiteva=zeros(ps,1);
minmse=min(posmse);


