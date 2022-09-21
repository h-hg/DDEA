% 'startSTR'
function  [newdata_x, newdata_f] = STR(FUN, Seed, CostValue)
    Range = [min(Seed); max(Seed)]';
    SeedSize = size(Seed,1);                     % Total Seedulation size
    D = size(Seed,2);
    SeedSize = size(Seed,1);
    dlimit=0.01;
    kernal = 'cubic';
    % Create the initial Seedulation
    CostValue = CostValue';
    [CostValue, inx] = sort(CostValue,'ascend');                            % Sort objective function in descending order
    Seed = Seed(inx,:);                                                     % Sort Initial Seedulation in descending order of its objective function
    global Dim
    Dim = D;
    
    %=====================================================================
    %% Begin the evolution loop
    %--------------------------------------
    sumSeed=Seed;
    YsumSeed=CostValue;
    Pbest=CostValue(1);
    newdata_x = [];
    newdata_f = [];
    
    %%%%%%%%%%%%%%evaluated points used for metamodeling
    sumMX=Seed;
    sumMY=CostValue;
    
    
    %%%%%%%%%%%%%%%%�?部搜�?
    index = unidrnd(1);
    Seedlocal=Seed(index,:);
    YSeedlocal=CostValue(index);
%     'startSTR'
    % Begin the evolution loop
    tr = 1;
    kmax = 3;
    for tr=1:kmax
        if tr==1
            %%%%%%%%%%%%%%%%%%%%%建立�?部代理模�? 
            DI=[];
            XI=[];
            YI=[];
            for j=1:length(sumMX(:,1))
                d=norm(Seedlocal-sumMX(j,:));
                if d>dlimit
                    dj=fdistance(sumMX(j,:),XI);
                    if dj>dlimit
                        XI=[XI;sumMX(j,:)] ; 
                        YI=[YI;sumMY(j,:)];
                        DI=[DI;d];
                    end
                end
            end
            ModelX=[];YModelX=[];
            len=5*D;%%%%%%%%%%%%%%%%%%%%%%%%%%%%the number of evaluated points used for local metamodeling
            if size(DI,1)<len
                ModelX=XI;
                YModelX=YI;
            else
                [DI,POS]=sort(DI);
                XI=XI(POS,:);
                YI=YI(POS,:);
                ModelX=XI(1:len,:);
                YModelX=YI(1:len,:);
            end
            ModelX=[ModelX;Seedlocal];
            YModelX=[YModelX;YSeedlocal]; 
            [~,xmin]=min(YModelX);
            Xmin=ModelX(xmin,:);
            [~,xmax]=max(YModelX);
            Xmax=ModelX(xmax,:);
            radius=norm(Xmin-Xmax)/2;
        else
            lb=(Seedlocal-radius)';%%%%%%%%%%%%%%%%%%%local optimizaiton range
            ub=(Seedlocal+radius)';
            for j=1:D
                if lb(j,1)<Range(j,1)
                    lb(j,1)=Range(j,1);
                end
                if ub(j,1)>Range(j,2)
                    ub(j,1)=Range(j,2);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%
            DI=[];
            XI=[];
            YI=[];
            for j=1:length(sumMX(:,1))
                d=norm(Seedlocal-sumMX(j,:));
                if d>dlimit
                    dj=fdistance(sumMX(j,:),XI);
                    if dj>dlimit
                        XI=[XI;sumMX(j,:)] ; 
                        YI=[YI;sumMY(j,:)] ; 
                        DI=[DI;d];
                    end
                end
            end
            ModelX=[];YModelX=[];
            len=5*D;%%%%%%%%%%%%%%%%%%%%%%%%
            if size(DI,1)<len
                ModelX=XI;
                YModelX=YI;
            else
                [DI,POS]=sort(DI);
                XI=XI(POS,:);
                YI=YI(POS,:);
                tt=0;
                for j=1:D
                    tt=tt+(XI(len,j)<ub(j)&&XI(len,j)>lb(j));
                end
                if tt<D
                    ModelX=XI(1:len,:);
                    YModelX=YI(1:len,:);
                else
                    for k=len:size(XI,1)
                        tt=0;
                        for j=1:D
                            tt=tt+(XI(k,j)<ub(j)&&XI(k,j)>lb(j));
                        end
                        if tt<D
                            break;
                        end
                    end
                    ModelX=XI(1:k,:);
                    YModelX=YI(1:k,:);
                end
            end
            ModelX=[ModelX;Seedlocal];
            YModelX=[YModelX;YSeedlocal]; 
        end
        
        lb=(Seedlocal-radius)';
        ub=(Seedlocal+radius)';
        for j=1:D
            if lb(j,1)<Range(j,1)
                lb(j,1)=Range(j,1);
            end
            if ub(j,1)>Range(j,2)
                ub(j,1)=Range(j,2);
            end
        end
        global coefC
        coefC=[];
        coefC=rbfcreate(ModelX',YModelX','RBFFunction', kernal);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
%             AlgStrategy = unidrnd(2);
%             if AlgStrategy == 1
%                 Alg = 'quasi-newton';  
%             elseif AlgStrategy == 2
%                 Alg = 'trust-region';
%             end

            AlgStrategy = 2;
            if AlgStrategy == 1
                Alg = 'quasi-newton';  
            elseif AlgStrategy == 2
                Alg = 'trust-region';  
            end
            opts = optimset('Algorithm',Alg,'Display', 'off','MaxIter',400,'MaxFunEvals',40*Dim);
            number=9;
            Seedinitial=lhsdesign(number,Dim);
            for j=1:Dim
               Seedinitial(:,j)=Seedinitial(:,j)*(ub(j)-lb(j))+lb(j); 
            end
            Xnumber=[];Fnumber=[];
            for j=1:number   
                [x1,f1]=fminunc(@YCRBF,Seedinitial(j,:),[],[],[],[],lb,ub,[],opts);
                Xnumber=[Xnumber;x1];Fnumber=[Fnumber;f1];
            end
            [fap,id]=min(Fnumber);
            Nap=Xnumber(id,:);
            [Nap1,fap1]=fminunc(@YCRBF,Seedlocal,[],[],[],[],lb,ub,[],opts);
            if fap1<fap
                Nap=Nap1;
                fap=fap1;
            end
        catch
            [Nap,fap]=YPSO(@YCRBF,[lb,ub]);
            fap = YCRBF(Nap);
        end
        
        candidate_position = Nap;
        
        %%%%%%%%%%%%%%%%%判断信赖域的更新准则
        if distance(Nap,ModelX)>10^-10
            fap=FUN(Nap);
            newdata_x = [newdata_x; Nap];
            newdata_f = [newdata_f; fap];
            ro=(YSeedlocal-fap)/(YCRBF(Seedlocal)-YCRBF(Nap)-10^(-20));
            if fap<YSeedlocal
                Seedlocal=Nap;
                YSeedlocal=fap;
            end
            sumSeed=[sumSeed;Nap];
            
            Pbest=[Pbest;min(min(Pbest),fap)];
            %%%%%%%%%建模点保�?
            if min(pdist2(Nap,sumMX))>dlimit
                sumMX=[sumMX;Nap];
                sumMY=[sumMY;fap];
            else
                [~,index]=min(pdist2(Nap,sumMX));
                if fap<min(sumMY)
                    sumMX(index,:)=Nap;
                    sumMY(index)=fap;
                end
            end
        else
            break;
        end
        
        if norm(Seedlocal-Nap)<radius
            es=1;
        end
        if norm(Seedlocal-Nap)>=radius
            es=2;
        end
        radius1 = radius;
        if ro<0.25
            radius1=0.25*radius;
        end
        if ro>0.75
            radius1=es*radius;
        end
        if ro>=0.25&&ro<=0.75
            radius1=radius;
        end
        radius=radius1;
    end
    
    if isempty(newdata_f)
        newdata_x = [newdata_x; Nap];
        newdata_f = [newdata_f; fap];
    end
end