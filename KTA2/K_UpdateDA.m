function [DAobj,DAdec,DAvar] = K_UpdateDA(DAobj,DAdec,DAvar,MaxSize,p)
% Update DA

    DA_Nor_pre = (DAobj - repmat(min(DAobj,[],1),size(DAobj,1),1))...
./repmat(max(DAobj,[],1) - min(DAobj,[],1),size(DAobj,1),1);

    %% Find the non-dominated solutions

    ND = NDSort(DAobj,1);
    DAobj = DAobj((ND==1),:);
    DAdec = DAdec((ND==1),:);
    DAvar = DAvar((ND==1),:);
    DA_Nor_pre = DA_Nor_pre((ND==1),:);
    N  = size(DAobj,1);
    if N <= MaxSize
        return;
    end
    
    %% Select the extreme solutions first
    
    Choose = false(1,N);
    
    M = size(DA_Nor_pre,2);
    select = randperm(M);
     
   
%     for i = 1:M
%    index =  find(associate == i);
%    [~,index1(i)] = min(val(index));
%     end
    
   % [~,Extreme1] = min(DA_Nor_pre,[],1);
    %[~,Extreme2] = max(DA_Nor_pre,[],1);
    Choose(select(1)) = true;
   % Choose(Extreme2) = true;
    
    
    
    
    

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
                Distance(i,j) = norm(DA_Nor_pre(i,:)-DA_Nor_pre(j,:),p);
                Distance(j,i) = Distance(i,j);
            end
        end
        while sum(Choose) < MaxSize
            Remain = find(~Choose);
            [~,x]  = max(min(Distance(~Choose,Choose),[],2));
            Choose(Remain(x)) = true;
        end
    end
%     DA = DA(Choose);
    DAobj = DAobj(Choose,:);
    DAdec = DAdec(Choose,:);
    DAvar = DAvar(Choose,:);
end
