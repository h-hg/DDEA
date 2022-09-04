function CA = UpdateCA(CA,New,MaxSize)
% Update CA
    if isempty(CA)
       CAObj = New.obj;
       CADec = New.dec; 
    else
    CAObj = [CA.obj;New.obj];
    CADec = [CA.dec;New.dec];
    end
    N = size(CAObj,1);
    
%     CA = [CA,New];
%     N  = length(CA);
    if N <= MaxSize
        CA.obj = CAObj;
        CA.dec = CADec;
        return;
    end
    
    %% Calculate the fitness of each solution
%     CAObj = CA.obj;
    CAObj = (CAObj-repmat(min(CAObj),N,1))./(repmat(max(CAObj)-min(CAObj),N,1));
    I = zeros(N);
    for i = 1 : N
        for j = 1 : N
            I(i,j) = max(CAObj(i,:)-CAObj(j,:));
        end
    end
    C = max(abs(I));
    F = sum(-exp(-I./repmat(C,N,1)/0.05)) + 1;
    
    %% Delete part of the solutions by their fitnesses
    Choose = 1 : N;
    while length(Choose) > MaxSize
        [~,x] = min(F(Choose));
        F = F + exp(-I(Choose(x),:)/C(Choose(x))/0.05);
        Choose(x) = [];
    end
    
    CA.obj = CAObj(Choose,:);
    CA.dec = CADec(Choose,:);
 
end