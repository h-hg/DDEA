function DA = UpdateDA(DA,New,MaxSize,p)
% Update DA

    %% Find the non-dominated solutions
    DAObj = [DA.obj;New.obj];
    DADec = [DA.dec;New.dec];
    
    ND = NDSort(DAObj,1);
    DAObj = DAObj(ND==1,:);
    DADec = DADec(ND==1,:);  
    N  = size(DAObj,1);
    if N <= MaxSize
        DA.obj = DAObj;
        DA.dec = DADec;
        return;
    end
    
    %% Select the extreme solutions first
    Choose = false(1,N);
    [~,Extreme1] = min(DAObj,[],1);
    [~,Extreme2] = max(DAObj,[],1);
    Choose(Extreme1) = true;
    Choose(Extreme2) = true;
    
    %% Delete or add solutions to make a total of K solutions be chosen by truncation
    if sum(Choose) > MaxSize
        % Randomly delete several solutions
        Choosed = find(Choose);
        k = randperm(sum(Choose),sum(Choose)-MaxSize);
        Choose(Choosed(k)) = false;
    elseif sum(Choose) < MaxSize
        % Add several solutions by truncation strategy
        Distance = inf(N);
        for i = 1 : N-1
            for j = i+1 : N
                Distance(i,j) = norm(DAObj(i,:)-DAObj(i,:),p);
                Distance(j,i) = Distance(i,j);
            end
        end
        while sum(Choose) < MaxSize
            Remain = find(~Choose);
            [~,x]  = max(min(Distance(~Choose,Choose),[],2));
            Choose(Remain(x)) = true;
        end
    end
    DA.obj = DAObj(Choose,:);
    DA.dec = DADec(Choose,:);
end