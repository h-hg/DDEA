function [ gbest,Record,time ] =CALSAPSO( Data,BU,BD,problem_name)
% Usage: [ gbest,Record,time ] =CALSAPSO( Data,BU,BD,problem_name)
% -----------------------------------------------------------------
% Important Note: This code needs intralled SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:
% Data          - Offline Data with c Decision Variables and Exact Objective Value
% problem_name  - Name of problem
% BU            - Upper Boundary of c Decision Variables
% BD            - Lower Boundary of c Decision Variables
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
warning off all;
tic;
c=size(Data,2)-1;
n=size(Data,1);
Record=[];
n1=0;
n2=0;
n3=0;

%--------StepControl:1V 2D 3L-------------------

StepNext=1;
Tm=min(Data(:,c+1));

while size(Data,1)<11*c
    StepCur=StepNext;
    switch StepCur
        case 1
            P=sortrows(Data,c+1);
            DataT= selectData( P,n+n1,c+1,2);
            DataT=sortrows(DataT,c+1);
            DataT(1,:)=P(1,:);
            [time,gbest,POP ] =PSO_DV(DataT,BU,BD);
        case 2
            P=sortrows(Data,c+1);
            DataT= selectData( P,n+n1,c,2);
            DataT=sortrows(DataT,c+1);
            DataT(1,:)=P(1,:);
            [time,gbest,POP ] =PSO_D(DataT,BU,BD);
        case 3
            P=sortrows(Data,c+1);      
            I=find(P(:,c+1)<Tm);
            DataT=P(I,1:c+1);
            if size(DataT,1)<3
                DataT= selectData( P,n+n1,c+1,2);
                DataT=sortrows(DataT,c+1);
                DataT(1,:)=P(1,:);
                [time,gbest,POP ] =PSO_DV(DataT,BU,BD);
            else
                bu=max(DataT(:,1:c),[],1);
                bd=min(DataT(:,1:c),[],1);
                [time,gbest,POP ] =PSO_DB(DataT,bu,bd);
            end
    end
    gbest(:,c+1)=compute_objectives(gbest(:,1:c),c,problem_name);
    t=[gbest,StepCur];
    Record=[Record;t];
    switch StepCur
        case 1
            StepNext=2;
            n3=n3+1;
        case 2
            if min(Data(:,c+1))>=min(gbest(:,c+1))
                StepNext=1;
                n1=n1+1;
            else
                StepNext=3;  
            end
        case 3
            if min(Data(:,c+1))>=min(gbest(:,c+1))
                StepNext=3;
                n2=n2+1;
            else
                StepNext=1;
            end
    end
    Data=[Data;gbest(:,1:c+1)];
    StepPre=StepCur;
end

[A,I]=min(Data(:,c+1));
gbest=Data(I,:);
toc;
time=toc;
end

