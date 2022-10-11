function [PopDec,PopObj] = Sur_MOEAD(PopDec,PopObj,Ensemble,r,Global,W)

W  = UniformPoint(45,Global.M);
Model1 = Ensemble{1};  Model2 = Ensemble{2};  Model3 = Ensemble{3};
r1 = r{1};  r2 = r{2};  r3 = r{3};
[f,~] = NDSort(PopObj,1);
fin = find(f==1);
rin = find(f~=1);
FObj = PopObj(fin,:);
FDec = PopDec(fin,:);
RDec = PopDec(rin,:);
RObj = PopObj(rin,:);
% if size(PopDec,1) > 45
%     if length(fin) < 45
%         [ind,~]= kmeans(RObj,45-length(fin));
%         Next   = zeros(1,45-length(fin));
%         for i = unique(ind)'
%             current = find(ind==i);
%             if length(current)>1
%                 best = randi(length(current),1);
%             else
%                 best = 1;
%             end
%             Next(i)  = current(best);
%         end
%         PopDec = [FDec;RDec(Next,:)];
%         PopObj = [FObj;RObj(Next,:)];
%     else
%         [ind,~]= kmeans(FObj,45);
%         Next   = zeros(1,45);
%         for i = unique(ind)'
%             current = find(ind==i);
%             if length(current)>1
%                 best = randi(length(current),1);
%             else
%                 best = 1;
%             end
%             Next(i)  = current(best);
%         end
%         PopDec = FDec(Next,:);
%         PopObj = FObj(Next,:);
%     end
% end
% ind = randperm(size(PopDec,1));
% PopDec = PopDec(ind,:);
% PopObj = PopObj(ind,:);

NDec = PopDec(fin,:);
NObj = FObj;
if length(fin) < 10
    NDec = [FDec;RDec(randperm(size(RDec,1),10-length(fin)),:)];
    NObj = [FObj;RObj(randperm(size(RDec,1),10-length(fin)),:)];
end
% if length(fin) > 45
%     [ind,~]= kmeans(FObj,45);
%         Next   = zeros(1,45);
%         for i = unique(ind)'
%             current = find(ind==i);
%             if length(current)>1
%                 best = randi(length(current),1);
%             else
%                 best = 1;
%             end
%             Next(i)  = current(best);
%         end
%         NDec = FDec(Next,:);
%         NObj = FObj(Next,:);
% end
ind = randperm(size(NDec,1));
PopDec = NDec(ind,:);
PopObj = NObj(ind,:);

W  = UniformPoint(size(PopObj,1),Global.M);
if size(W,1)>size(PopObj,1)
    W = W(randperm(size(W,1),size(PopObj,1)),:);
end
size(W,1)
g = 0;  gmax = 20;  N = size(PopDec,1);
T = ceil(N/10);
if T < 2
    T = 2;
end
B = pdist2(W,W);
[~,B] = sort(B,2);
B = B(:,1:T);
Z = min(PopObj,[],1);

while g <= gmax
    for k = 1 : size(W,1)
        P = B(k,randperm(size(B,2)));
        OffDec = GAhalf(PopDec(P(1:2),:));
        for i = 1: size(OffDec,1)
            for j = 1 : Global.M
                [y1(i,j),~,MSE1(i,j)] = predictor(OffDec(i,r1),Model1{j});
                [y2(i,j),~,MSE2(i,j)] = predictor(OffDec(i,r2),Model2{j});
                [y3(i,j),~,MSE3(i,j)] = predictor(OffDec(i,r3),Model3{j});
            end
        end
        theta = Global.evaluated/Global.evaluation;
%         OffObj = 0.5*(1-theta)*(y1+y2)/2 + 0.5*(1+theta)*y3;
        OffObj = theta * y3 + (1 - theta) * (y1+y2)/2;
        Z = min(Z,OffObj);
        normW   = sqrt(sum(W(P,:).^2,2));
        normP   = sqrt(sum((PopObj(P,:)-repmat(Z,T,1)).^2,2));
        normO   = sqrt(sum((OffObj-Z).^2,2));
        CosineP = sum((PopObj(P,:)-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
        CosineO = sum(repmat(OffObj-Z,T,1).*W(P,:),2)./normW./normO;
        g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
        g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
        nr = sum(g_old>=g_new);
        PopDec(P(g_old>=g_new),:) = repmat(OffDec,nr,1);
        PopObj(P(g_old>=g_new),:) = repmat(OffObj,nr,1);
    end
    g = g+1;
end
end

