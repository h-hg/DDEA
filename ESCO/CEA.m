function [Best,Data0] = CEA(Data,BU,BD,problem_name)

dim = size(Data,2)-1;
num0 = size(Data,1);
pc=1.0;  %Crossover Probability 交叉概率
pm=1/dim;  %Mutation Probability 变异概率
bu = BU;  bd = BD;
% T = 100;
% K = ceil(T*0.1);
T = 10;
K = ceil(T*0.3);  
l = sqrt(0.001^2*dim);
N = 50;
gmax = 20; %100;
FEs = 1000;
fes = num0;
Data0 = Data;

while fes <= FEs
    fprintf('当前函数计算: fes = %d 次\n', fes);
    Data = Data_Process(Data0,num0);
    num = size(Data,1);
    Best = min(Data(:,dim+1));
    Cr = KRCC(Data,num,dim);
    [Model,Sind,Xind,Train,Test,W,B,C,P] = Low_dim_RBF(Data,num,dim,Cr,T);
    [Ens,Mind,Sim] = Selective_Ensemble(Model,Data,Train,Test,K,W,B,C,P);
    
    [X,Y] = Sur_Coevolution(Data,Ens,Train,Xind,Mind,Sim,bu,bd,N,gmax,W,B,C,P);
    new = Infill_solution_Selection(Data,[X,Y],l);
    
    if ~isempty(new)
        newY = compute_objective(new(:,1:end-1),dim,problem_name);
        Data0 = [ Data0;[new(:,1:end-1),newY] ];
        fes = fes + length(newY);
    end
    
end  %%end while
Best = min(Data0(:,dim+1));

end  %%end function

function Data = Data_Process(Data,num0)
num = size(Data,1);
if num > num0+400
    [~,ind] = sort(Data(:,end));
    c = num0+400+1;
    Data(ind(c:end),:) = [];
end
end
function Cr = KRCC(S,n,d)

X = S(:,1:d);
Y = S(:,d+1);
corr = zeros(1,d);

for i = 1 : d
    nc = 0; nd = 0;
    for j = 1 : n
        for k = j+1 : n
            if X(j,i) > X(k,i) && Y(j,1) > Y(k,1) || X(j,i) < X(k,i) && Y(j,1) < Y(k,1)
                nc = nc + 1;
            else
                nd = nd + 1;
            end
        end
    end
    corr(i) = (nc - nd)/( n*(n-1)/2 );
end
Cr = abs(corr)/sum(abs(corr));
end
function y = Ens_predictor(X,Ens,Train,Xind,Mind,Sim)

N = size(X,1);
K = length(Ens);
pre = zeros(N,K);
for i = 1 : K
    
    S = Train{Mind(i)};
    pre(:,i) = rbfpredict(Ens{i},S(:,1:end-1),X(:,Xind{Mind(i)}));
end
% y = mean(pre,2);
sim = 1./Sim;
w = sim(Mind(1:K))/sum(sim(Mind(1:K)));
wpre = repmat(w,N,1) .* pre;
y = sum(wpre,2);
end

% rng('shuffle');
% k = randperm(length(Ens),1);
% model = Ens{k};  xind = Xind{Mind(k)};
% ldim = length(xind);  %辅助子搜索的维数
% train = Train{Mind(k)};
% X = initialize_pop(N,dim,bu,bd); %初始化辅助搜索的初始种群aX
% Y = Ens_predictor(X,Ens,Train,Xind,Mind,Sim);
% aX = initialize_pop(N,ldim,bu(:,xind),bd(:,xind)); %初始化辅助搜索的初始种群aX
% aY = rbfpredict(model,train(:,1:end-1),aX);
% g = 0;
% while g <= 500
%     [Ens,Mind,Sim] = Selective_Ensemble(Model,Data0,Train,Test,K);  
%     X1 = SBX(X,bu,bd,pc,N/2);
%     X2 = mutation(X1(1:N/2,:),bu,bd,pm,N/2);
%     aX1 = SBX(aX,bu(:,xind),bd(:,xind),pc,N/2);
%     aX2 = mutation(aX1(1:N/2,:),bu(:,xind),bd(:,xind),pm,N/2);
%     X3 = X2;  X3(:,xind) = aX2;
%     X = [X;X2;X3];
%     aX3 = X2(:,xind);
%     aX = [aX;aX2;aX3];
%     y = Ens_predictor([X2;X3],Ens,Train,Xind,Mind,Sim);
%     ay = rbfpredict(model,train(:,1:end-1),[aX2;aX3]);
%     Y = [Y;y];
%     aY = [aY;ay];
%     %     y = Ens_predictor([X2;X3],Ens,Train,Xind,Mind);
%     %     Y = [Y;y];
%     %     Y = Ens_predictor(X,Ens,Train,Xind,Mind);
%     %     aY = rbfpredict(model,train(:,1:end-1),aX);
%     [sY,is1] = sort(Y);  %适应值排序
%     [saY,is2] = sort(aY);
%     X = X(is1(1:N),:);  %选择前N个
%     Y = Y(is1(1:N),1);
%     aX = aX(is2(1:N),:);
%     aY = aY(is2(1:N),1);
%     
%     new = X(1,:);
%     obj = compute_objective(new,dim,problem_name);
%     Data = [Data;[new,obj]];
%     g = g + 1;
% end

