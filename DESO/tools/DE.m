function [bestP,bestFitness] = DE(Dim, Max_NFEs,FUN,minerror,ghx)
    time_begin=tic;
    n = Dim; NP = 50; flag_er=0;
    lu = [min(ghx); max(ghx)];
    LowerBound = lu(1, :);
    UpperBound = lu(2, :);
    
    %% DE Parameters
    G=1;
    F=0.5; 
    CR=0.9;
    UB=UpperBound;
    LB=LowerBound;
    
    UB=repmat((UB),NP,1);
    LB=repmat((LB),NP,1);
    
    P=(UB-LB).*rand(NP,Dim)+LB;
    
    fitnessP=FUN(P);
        
    NFEs=NP;
    [fitnessBestP,indexBestP]=min(fitnessP);
    bestP=P(indexBestP,:);
    recRMSE(1:NP)=fitnessP;
    
    while NFEs<Max_NFEs
        
        fitnessBestP_old = fitnessBestP;
        for i=1:NP  
            k0=randi([1,NP]);
            while(k0==i)
                k0=randi([1,NP]);   
            end
            P1=P(k0,:);
            k1=randi([1,NP]);
            while(k1==i||k1==k0)
                k1=randi([1,NP]);
            end
            P2=P(k1,:);
            k2=randi([1,NP]);
            while(k2==i||k2==k1||k2==k0)
                k2=randi([1,NP]);
            end
            P3=P(k2,:);
            
            V(i,:)=P1+F.*(P2-P3);   
            
            for j=1:Dim
              if (V(i,j)>UB(i,j)||V(i,j)<LB(i,j))
                 V(i,j)=LB(i,j)+rand*(UB(i,j)-LB(i,j));         
              end
            end
            
            jrand=randi([1,Dim]); 
            for j=1:Dim
                k3=rand;
                if(k3<=CR||j==jrand)
                    U(i,j)=V(i,j);
                else
                    U(i,j)=P(i,j);
                end
            end
        end
        
        time_begin=tic;
        fitnessU=FUN(U);
        time_cost=toc(time_begin);
        
        NFEs=NFEs+NP;
        
        for i=1:NP 
            if(fitnessU(i)<fitnessP(i))
                P(i,:)=U(i,:);
                fitnessP(i)=fitnessU(i);
                if(fitnessU(i)<fitnessBestP)
                   fitnessBestP=fitnessU(i);
                   bestP=U(i,:);
                end
            end
            recRMSE(NFEs)=fitnessP(i);
        end
        error=abs(fitnessBestP_old-fitnessBestP); 
        if error <= minerror   
            flag_er=flag_er+1;
        else
            flag_er=0;
        end
        if flag_er >=10
            break;
        end
        G=G+1;
    end
    bestFitness=fitnessBestP;
    endNFEs = NFEs;
    time_cost=toc(time_begin);
end
    