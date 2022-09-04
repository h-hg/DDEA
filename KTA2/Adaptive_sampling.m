function Offspring01 = Adaptive_sampling(CAobj,DAobj,CAdec,DAdec,DAvar,DA,mu,p,phi)

Ideal_Point = min([CAobj;DAobj],[],1);  

flag = Cal_Convergence(CAobj,DAobj,Ideal_Point); 

if flag == 1   
    % convergence sampling strategy
    N = size(CAobj,1);
    CAobj01 = (CAobj-repmat(min(CAobj),N,1))./(repmat(max(CAobj)-min(CAobj),N,1));
    I = zeros(N);
    for i = 1:N
        for j = 1:N
            I(i,j) = max(CAobj01(i,:)-CAobj01(j,:));
        end
    end
    C = max(abs(I));
    F = sum(-exp(-I./repmat(C,N,1)/0.05)) + 1;
    Choose = 1:N;
    while length(Choose) > mu
        [~,x] = min(F(Choose));
        F = F + exp(-I(Choose(x),:)/C(Choose(x))/0.05);
        Choose(x) = [];
    end
    Offspring01 = CAdec(Choose,:);
else
    if PD(DAobj,[]) < PD(DA.obj,[])
        % uncertainty sampling strategy
        An = size(DAvar,1);
        Choose = zeros(1,5);
        for i = 1:mu
            A_num = randperm(size(DAvar,1));
            Uncertainty = mean(DAvar(A_num(1:ceil(phi*An)),:),2);
            [~,best]    = max(Uncertainty);
            Choose (i)     = A_num(best);
        end
        Offspring01 = DAdec(Choose ,:);
    else
        % diversity sampling strategy
        DA_Nor = (DA.obj - repmat(min([DAobj;DA.obj],[],1),length(DA),1))...
            ./repmat(max([DAobj;DA.obj],[],1) - min([DAobj;DA.obj],[],1),size(DA.obj,1),1);
        DA_Nor_pre = (DAobj - repmat(min([DAobj;DA.obj],[],1),size(DAobj,1),1))...
            ./repmat(max([DAobj;DA.obj],[],1) - min([DAobj;DA.obj],[],1),size(DAobj,1),1);
        N  = size(DA_Nor,1);
        Pop = [DA_Nor;DA_Nor_pre];
        Pop_dec = [DA.dec;DAdec];
        NN = size(Pop,1);
        Choose = false(1,NN);
        Choose(1:N) = true;
        MaxSize = N+mu;
        Distance = inf(N);
        for i = 1 : NN-1
            for j = i+1 : NN
                Distance(i,j) = norm(Pop(i,:)-Pop(j,:),p);
                Distance(j,i) = Distance(i,j);
            end
        end
        Offspring01 = [];
        while sum(Choose) < MaxSize
            Remain = find(~Choose);
            [~,x]  = max(min(Distance(~Choose,Choose),[],2));
            Choose(Remain(x)) = true;
            Offspring01 = [Offspring01;Pop_dec(Remain(x),:)];
        end
    end
end



end