function [ time,gbest,logbook] =TSABFEA( problem_name,c)
% Usage: [ time,gbest,logbook] =TSABFEA( problem_name,c)
%
% Input:
% problem_name  - Benchmark 
% c             - No. of Decision Variables
%
% Output: 
% time          - Execution Time
% gbest         - Obtained Optimum
% logbook       - Search history
%
    %%%%    Authors:    Handing Wang, Yaochu Jin,Cuie Yang, Licheng Jiao
    %%%%    Xidian University, China, University of Surrey, UK and Northeastern University, China
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage/home
    %%%%    DATE:       August 2020
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:

%Handing Wang, Yaochu Jin, Cuie Yang, and Licheng Jiao, Transfer stacking from low-to high-fidelity: A surrogate-assisted bi-fidelity evolutionary algorithm, Applied Soft Computing, accepted.

%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
rand('state',sum(100*clock));

switch problem_name
    case {'MFB1','MFB2','MFB3','MFB7','MFB8','MFB9','MFB10','MFB11','MFB12','MFB13'}
        type='C';
        Bphi=[1000,10000];
    case {'MFB4'} 
        type='D';
        Bphi=[0:1000:10000];
    case {'MFB5'} 
        type='D';
        Bphi=[1000,3000,10000];
    case {'MFB6'} 
        type='D';
        Bphi=[1000,10000];
    otherwise
        type='C';
        Bphi=[1000,10000];
end
PHImax=max(Bphi);
phi=min(Bphi);
bu=ones(1,c);
bd=0-ones(1,c);


logbook=[];
Cost=0;
n=100;
nc=c;
tic;
TrainSizemax=n;
% TrainSizemin=0.5*n;
TrainSize=TrainSizemax;
%=============Train RBF==============
POP = initialize_popLHS(n,c,bu,bd);
[ obj,f,e,cost ] =compute_objectives(POP,c,phi,problem_name);
Cost=Cost+sum(cost);
POP=[POP,obj];
[ obj,f,e,cost ] =compute_objectives(POP,c,PHImax,problem_name);
Cost=Cost+sum(cost);
Train=[POP,obj];

%=============POP initialization==============
POP = initialize_pop(n,c,bu,bd);
[ obj,f,e,cost ] =compute_objectives(POP,c,phi,problem_name);
Cost=Cost+sum(cost);
POP=[POP,obj];


pc=1;
pm=1/c;
g=0;
gmax=20000;
Costmax=10000000;
nls=2;
[A,Imin]=min(Train(:,c+2));
gbest=Train(Imin,1:c+2);
t=[Cost,gbest];
logbook=[logbook;t];
while  Cost<Costmax&g<gmax
    NPOP= SBX( POP,bu,bd,pc,n );
    NPOP=mutation(NPOP,bu,bd,pm);

    [ obj,f,e,cost ] =compute_objectives(NPOP,c,phi,problem_name);
    NPOP=[NPOP,obj];
    Cost=Cost+sum(cost);
    POP=[POP;NPOP];
    [ W2,B2,Centers,Spreads] = RBF( Train(:,1:c),Train(:,end),nc);
    obj=RBF_predictor(W2,B2,Centers,Spreads, POP(:,1:c));
    POP=[POP,obj];

    IN= NearestNeighbor(POP,Train,c );
    obj=RBF_predictor(W2,B2,Centers,Spreads,Train(IN,1:c));
    
    A=[Train(IN,c+1),obj];
    w=A\Train(IN,end);
    obj=[POP(:,c+1:c+2)]*w;
    
    POP=[POP,obj];  

    [ obj,f,e,cost ] =compute_objectives(POP,c,phi,problem_name);
    POP=[POP,f];

    POP=sortrows(POP,c+3);
    if size(Train,1)>TrainSize
        Train=Train(size(Train,1)-TrainSize+1:end,:);
    end
   

    EX=POP;
    dis= MinDis(EX,Train,c );
    EX=EX(find(dis~=0),:);
    if isempty(EX)
        dis = MinDis(POP,Train,c );
        I=find(dis~=0);
        [A,Ii]=sort(0-abs(POP(I,c+1)*sign(w(1))-POP(I,c+2)*sign(w(2))));
        Ib=[I(Ii(1));I(1)];
        Ib=unique(Ib);
        EX=POP(Ib,:);
    elseif size(EX,1)>nls
        [A,Ii]=sort(0-abs(EX(:,c+1)-EX(:,c+2)));
        Ib=[1;Ii(1)];
        Ib=unique(Ib);
        EX=EX(Ib,:);
    end

    [ obj,f,e,cost ] =compute_objectives(EX,c,PHImax,problem_name);
    Cost=Cost+sum(cost);
    T=[EX(:,1:c+1),obj];
    Train=[Train;T];
    [A,Imin]=min(obj);
    gbest=EX(Imin,1:c+1);
    t=[Cost,gbest,f(Imin)];
    logbook=[logbook;t];
    
    T=f(Imin);

    [A,Imin]=min(logbook(:,end));
    gbest=logbook(Imin,2:c+3);
    if MinDis(gbest,POP(1:n,:),c )~=0
        dis=MinDis(Train,POP(1:n,:),c );
        I=find(dis~=0&Train(:,end)<T);
        POP=[gbest(1:c+1);Train(I,1:c+1);POP(:,1:c+1)];
    end
    POP=POP(1:n,1:c+1);
    g=g+1;
end
time=toc;
logbook(end,end);
end