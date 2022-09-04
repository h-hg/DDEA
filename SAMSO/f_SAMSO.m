%% ********************************************************************************************************************************************************8
% function  [fy,fNFE]=f_SAMSO(fun,D,lb,ub)
% Usage:  [fy,fNFE]=f_SAMSO(fun,D,lb,ub)
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
% Input:
% fun           - Name of the problem
% D              - Dimension of the problem
% lb             - Lower Boundary of Decision Variables
% ub            - Upper Boundary of Decision Variables
%
% Output:
% fy             - the obtained best fitness value
% fNFE         -  the obtained best fitness value in each iteration
%--------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------
% Authors:      Fan Li, Xiwen Cai, Liang Gao, Weiming Shen
% Address       Huazhong University of Science & Technology, Wuhan, PR China;
% WEBSITE:    https://sites.google.com/site/handingwanghomepage
% DATE:         April 2020
%This code is part of the program that produces the results in the following paper:
%Fan Li, Xiwen Cai, Liang Gao, Weiming Shen. A Surrogate-Assisted Multiswarm Optimization Algorithm for High-Dimensional Computationally Expensive Problems , IEEE Transactions on Cybernetics, 2020.
%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%% *********************************************************************************************************************************************************************************************************


function [fy,fNFE]=f_SAMSO(fun,D,lb,ub)
% lb=-20;ub=20;D=30; fun=@Rosenbrock;

Dat.myFN =fun;f=Dat.myFN;
Dat.designspace = [lb;ub]*ones(1,D);  % boundary
Dat.ndv = D;
bound=Dat.designspace;
vbound(2,:)=0.5*(bound(2,:)-bound(1,:));vbound(1,:)=-vbound(2,:);
maxNFE=1000;

%% create DOE
if D<=50
    k=40;%Sample number
else
    k=2*D;
end
Dat.npoints =k;Doe=@DOELHS ;
[Xtrain Ytrain]=my_initial(Dat,Doe);
NFE=k;
%% initial population
[val ind ]=sort(Ytrain);
if D<=50
    N=k;
    pop=Xtrain;fpop=Ytrain;pbest=pop;   fpbest=fpop;
else
    N=80;
    pop=Xtrain(ind(1:N),:);fpop=Ytrain(ind(1:N));pbest=pop;   fpbest=fpop;
end
gbest=Xtrain(ind(1),:);fgbest=fpbest(1);

[u,~]=my_initial(Dat,Doe);
u=0.25*u; v=u;
w=0.729;   c1=1.491;  c2=1.491;
%%
fNFE(1:NFE)=fgbest;Samp=[Xtrain];YS=[Ytrain];
dlta=min(sqrt(0.001^2*D),0.00005*sqrt(Dat.ndv)*min(Dat.designspace(2,:)-Dat.designspace(1,:)));
PAU=[];FPAU=[];%局部不确定点
tept=0;tepg=0;
t=1;
%% 开始迭代
while NFE< maxNFE
    G(t,:)=gbest;  FG(t)=fgbest;
    
    srgtOPT=srgtsRBFSetOptions(Samp,YS, @my_rbfbuild, [],'CUB', 0.0002,1);
    srgtSRGT = srgtsRBFFit(srgtOPT);
    L2 =@(x)my_rbfpredict(srgtSRGT.RBF_Model, srgtSRGT.P, x);
    FE=3000;
    options = optimset('Algorithm','interior-point','Display','off','MaxFunEvals',FE,'TolFun',1e-8,'GradObj','off'); % run interior-point algorithm
    L=Dat.designspace(1,:);U=Dat.designspace(2,:);
    x= fmincon(L2,gbest,[],[],[],[],L,U,[],options);
    dx=min(sqrt(sum((repmat(x,size(Samp,1),1)-Samp).^2,2)));
    if dx>dlta
        fx=Dat.myFN(x);    Samp=[ Samp;x]; YS=[YS;fx];
        NFE=NFE+1;  fNFE(NFE)=min(YS);
        if fx<fgbest
            fgbest= fx; gbest=x;PAU=[PAU;x];FPAU=[FPAU;fx];
            tepg=1;                        %记录全局最优信息
        else
            tepg=0;
        end
    end
    
    wn=1;
    
    if maxNFE<1000
        N2=2+round(N*((maxNFE-NFE)/maxNFE).^wn);N1=N-N2;
    else
        if NFE<1000
            maxNFE1=1000;
            N2=2+round(N*((maxNFE1-NFE)/maxNFE1).^wn);N1=N-N2;
        else
            N2=2;N1=N-N2;
        end
        
    end
    
    SEL=min(round(2*D),size(Samp,1));
    [nouse, seq] = sort(YS);
    Xsel = Samp(seq(1:SEL), :);
    xmean = mean(Xsel);
    % covariance matrix calculation
    C =  1/(SEL-1)*(Xsel - xmean(ones(SEL,1), :))'*(Xsel - xmean(ones(SEL,1), :));
    C = triu(C) + transpose(triu(C,1)); % enforce symmetry
    [B,D1] = eig(C);
    % limit condition of C to 1e20 + 1
    if max(diag(D1)) > 1e20*min(diag(D1))
        tmp = max(diag(D1))/1e20 - min(diag(D1));
        C = C + tmp*eye(D);
        [B, D1] = eig(C);         %R的每一列为特征向量
    end
    
    %% 构造引导的粒子
    
    for i=1:N
        if i<N1
            if rand<0.5
                v(i,:)=w*v(i,:)+((gbest-pop(i,:))*B).*rand(1,D)*c1*B'+((pbest(i,:)-pop(i,:))*B).*rand(1,D)*c2*B';%
            else
                v(i,:)=w*v(i,:)+rand(1,D)*c1.*(gbest-pop(i,:))+rand(1,D)*c2.*(pbest(i,:)-pop(i,:));%
            end
            for j=1:D
                if  v(i,j)<vbound(1,j)
                    v(i,j)=vbound(1,j);tept=tept+1;
                end
                if  v(i,j)>vbound(2,j)
                    v(i,j)=vbound(2,j);tept=tept+1;
                end
            end
            pop1(i,:)=pop(i,:)+v(i,:);
            for j=1:D
                if pop1(i,j)<bound(1,j)
                    pop1(i,j)=bound(1,j);
                end
                if  pop1(i,j)>bound(2,j)
                    pop1(i,j)=bound(2,j);
                end
            end
            newpop(i,:)=pop1(i,:);
            pop(i,:)=pop1(i,:);
        else
            A = 1:N;A(i)=[];  j = A(randi(N-1));
            Step = pop(i,:) - pop(j,:);
            if fpop(j)< fpop(i)
                Step = -Step;
            end
            newpop(i,:)= pop(i,:) + rand(1,D).*Step;
            
            for j=1:D
                if newpop(i,j)<bound(1,j)
                    newpop(i,j)=bound(1,j);
                end
                if  newpop(i,j)>bound(2,j)
                    newpop(i,j)=bound(2,j);
                end
            end
        end
    end
    
    %% 更新个体、全局最优
    srgtOPT=srgtsRBFSetOptions(Samp,YS, @my_rbfbuild, [], 'CUB',0.0002,1);
    srgtSRGT = srgtsRBFFit(srgtOPT);
    predy= my_rbfpredict(srgtSRGT.RBF_Model, srgtSRGT.P, newpop);
    
    tep=0;
    for i=1:N
        d=min(sqrt(sum((repmat(newpop(i,:),size(Samp,1),1)-Samp).^2,2)));
        if  predy(i)<fpop(i) & d>dlta
            tep= tep+1;
            fnewpop(i)=Dat.myFN( newpop(i,:));
            Samp=[ Samp; newpop(i,:)]; YS=[YS;fnewpop(i)];
            NFE=NFE+1; fNFE(NFE)=min(YS);
            if fnewpop(i)<fpop(i)
                pop(i,:) = newpop(i,:);fpop(i)=fnewpop(i);
            end
            if fpop(i)<fpbest(i)
                pbest(i,:)=pop(i,:);fpbest(i)= fpop(i);
            end
            if fpbest(i)<fgbest
                gbest= pbest(i,:);fgbest=fpbest(i);
            end
        end
    end
    
    if tep==0
        [val,ind]=min(predy);
        fnew=Dat.myFN( newpop(ind,:));
        Samp=[ Samp; newpop(ind,:)]; YS=[YS;fnew];
        NFE=NFE+1; fNFE(NFE)=min(YS);
        if fnew<fpop(ind)
            pop(ind,:) = newpop(ind,:);fpop(ind)=fnew;
            if fpop(ind)<fpbest(ind)
                pbest(ind,:) = pop(ind,:);fpbest(ind)=fpop(ind);
            end
            if fnew<fgbest
                gbest= newpop(ind,:);fgbest=fnew;
            end
        end
    end
    
    t=t+1;
    NFE
end
fy=fgbest;
