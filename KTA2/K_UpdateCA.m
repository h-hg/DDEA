function [CAobj,CAdec,CAvar] = K_UpdateCA(CAobj,CAdec,CAvar,MaxSize)
% Update CA

    N  = size(CAobj,1);
    if N <= MaxSize
        return;
    end
    
    %% Calculate the fitness of each solution
%     CAObj = CA.objs;


    CAobj1 = (CAobj-repmat(min(CAobj),N,1))./(repmat(max(CAobj)-min(CAobj),N,1));
    I = zeros(N);
    for i = 1 : N
        for j = 1 : N
            I(i,j) = max(CAobj1(i,:)-CAobj1(j,:));
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
    CAobj = CAobj(Choose,:);
    CAdec = CAdec(Choose,:);
    CAvar = CAvar(Choose,:);
end
