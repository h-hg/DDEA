function [bestP,bestFitness] = JADE(Dim, Max_NFEs,FUN,minerror,ghx)
    addpath('./CauchyFunction');
    
    time_begin=tic;
%     disp('JADE global search');
    n = Dim; NP = 100; flag_er=0;
    lu = [min(ghx); max(ghx)];
    LowerBound = lu(1, :);
    UpperBound = lu(2, :);

    G=1;
    uCR=0.5;
    uF=0.5;
    p0=0.05;
    top=round(p0*NP);
    A=[];
    t=1;
    c=0.1;
    
    UB=UpperBound;
    LB=LowerBound;
    
    UB=repmat((UB),NP,1);
    LB=repmat((LB),NP,1);
    
    P=(UB-LB).*rand(NP,Dim)+LB;%随机产生初始种群个体
    
    fitnessP=FUN(P);%计算种群个体适应值
    
    NFEs=NP;%记录适应度函数调用次数
    [fitnessBestP,indexBestP]=min(fitnessP);
    bestP=P(indexBestP,:);
    recRMSE(1:NP)=fitnessP;
    
    %% 总体大循环
    while NFEs<Max_NFEs
        fitnessBestP_old = fitnessBestP;
        
        Scr=[];%初始成功参加变异的个体的交叉概率为空集
        Sf=[];%初始成功参加变异的个体的缩放因子为空集
        n0=1;%记录Scr，Sf中的元素个数
        
        %根据个体适应值进行排序，得到最好的前top个个体
        [~,indexSortP]=sort(fitnessP);
        
        for i=1:NP
            
            %对第G代的每个个体计算交叉概率和缩放因子
            CR(i)=normrnd(uCR,0.1);%服从正态分布
            F(i)=cauchyrnd(uF,0.1);%服从柯西分布

            if(CR(i)>1)
               CR(i)=1;
            elseif CR(i)<0
               CR(i)=0; 
            end
            while (F(i)<=0)
                F(i)=cauchyrnd(uF,0.1);
            end
            if (F(i)>1)
                F(i)=1;
            end
            
        end
        
        %%变异操作+交叉操作
        for i=1:NP  
                                                                                                                                                                                                                                                                                                                                                                                  
            %%变异操作
            for j=1:top
               bestTopP(j,:)=P(indexSortP(j),:); 
            end
            %从top个个体中随机选出一个作为Xpbest
            k0=randperm(top,1);
            Xpbest=bestTopP(k0,:);
            %从当前种群P中随机选出P1
            k1=randi([1,NP]);
            P1=P(k1,:);
            while (k1==i)
                k1=randi([1,NP]);
                P1=P(k1,:); 
            end
            %从P∪A中随机选出P2
            PandA=[P;A];
            [num,~]=size(PandA);
            k2=randi([1,num]);
            P2=PandA(k2,:);
            while (k2==i||k2==k1)
                k2=randi([1,num]);
                P2=PandA(k2,:); 
            end
            V(i,:)=P(i,:)+F(i).*(Xpbest-P(i,:))+F(i).*(P1-P2);   
       
             %%交叉操作
            jrand=randi([1,Dim]); 
            for j=1:Dim
                k3=rand;
                if(k3<=CR(i)||j==jrand)
                    U(i,j)=V(i,j);
                else
                    U(i,j)=P(i,j);
                end
            end
            
        end
        
        
        %% 边界处理
        for i=1:NP
           for j=1:Dim
              if (U(i,j)>UB(i,j)||U(i,j)<LB(i,j))
                 U(i,j)=(UB(i,j)-LB(i,j))*rand+LB(i,j); 
              end
           end
        end
        
        fitnessU=FUN(U);
        NFEs=NFEs+NP;
        
        for i=1:NP
        %% 选择操作

            if(fitnessU(i)<fitnessP(i))
                A(t,:)=P(i,:);
                P(i,:)=U(i,:);
                fitnessP(i)=fitnessU(i);
                Scr(n0)=CR(i);
                Sf(n0)=F(i);
                t=t+1;
                n0=n0+1;
                if(fitnessU(i)<fitnessBestP)
                   fitnessBestP=fitnessU(i);
                   bestP=U(i,:);
                end
            end
            
            recRMSE(NFEs)=fitnessP(i);
        end

        %判断归档种群A的规模是否在NP之内，若大于，则随机移除个体使其规模保持NP
        [tA,~]=size(A);
        if tA>NP
            nRem=tA-NP;
            k4=randperm(tA,nRem);
            A(k4,:)=[]; 
            [tA,~]=size(A);
            t=tA+1;
        end

        %自适应参数，更新uCR和uF
        [~,ab]=size(Scr);
        if ab~=0
            newSf=(sum(Sf.^2))/(sum(Sf));
            uCR=(1-c)*uCR+c.*mean(Scr);
            uF=(1-c)*uF+c.*newSf;
        end
        %% 判别终止条件
        error=abs(fitnessBestP_old-fitnessBestP); 
        if error <= minerror   % minerror终止条件
            flag_er=flag_er+1;
        else
            flag_er=0;
        end
        if flag_er >=10
            break;
        end
        
%         disp(['Iteration ' num2str(G) '   Best Cost = ' num2str(fitnessBestP) '   NFE=' num2str(NFEs)]);
        G=G+1;
        
    end
    bestFitness=fitnessBestP;
    endNFEs = NFEs;
    
    time_cost=toc(time_begin);
end