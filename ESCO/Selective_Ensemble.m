function [Ens,Mind,Sim] = Selective_Ensemble(Model,Data,Train,Test,K,W,B,C,P)

Sim = zeros(1,length(Model));
Ens = cell(1,K);
for i = 1 : length(Model)
    
    S = Train{i};  Q = Test{i};
%     pre = rbfpredict(Model{i},S(:,1:end-1),Q(:,1:end-1));
    pre = RBF_predictor(W(i,:),B(i),C{i},P(:,i),Q(:,1:end-1));
    mpre = mean(pre);
    vpre = var(pre);
    my = mean(Data(:,end));
    vy = var(Data(:,end));
    sim = KLD(my,vy,mpre,vpre);
    Sim(i) = sim;  
end

[~,Mind] = sort(Sim);
Ens = Model(Mind(1:K));
% for i = 1 : K    
%     Ens{i} = Model(Mind(i));    
% end

end

function sim = KLD(mu1,sigma1,mu2,sigma2)

sim = log( sqrt(sigma2)/sqrt(sigma1) )+( sigma1+(mu1-mu2)^2 )/(2*sigma2)-1/2; 
end