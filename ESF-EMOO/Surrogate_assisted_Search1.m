function [PopDec,PopObj] = Surrogate_assisted_Search1(PopDec,PopObj,Ensemble,V0,r,Global)

V = V0;
Model1 = Ensemble{1};  Model2 = Ensemble{2};  Model3 = Ensemble{3};
r1 = r{1};  r2 = r{2};  r3 = r{3}; 
[f,~] = NDSort(PopObj,1);
fin = find(f==1);
rin = find(f~=1);
FObj = PopObj(fin,:);
FDec = PopDec(fin,:);
RDec = PopDec(rin,:);
RObj = PopObj(rin,:);
% if size(PopDec,1) > 50
%     if length(fin) < 50
%         [ind,~]= kmeans(RObj,50-length(fin));
%         Next   = zeros(1,50-length(fin));
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
%         [ind,~]= kmeans(FObj,50);
%         Next   = zeros(1,50);
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

NDec = PopDec(fin,:);
NObj = FObj;
if length(fin) < 10
    NDec = [FDec;RDec(randperm(size(RDec,1),10-length(fin)),:)];
    NObj = [FObj;RObj(randperm(size(RDec,1),10-length(fin)),:)];
end
ind = randperm(size(NDec,1));
PopDec = NDec(ind,:);
PopObj = NObj(ind,:);  

g = 0;  gmax = 20;
Off2 = GA(PopDec);
Off = [PopDec;Off2];
while g <= gmax
    drawnow();
    N     = size(Off2,1);
    y1 = zeros(N,Global.M); y2 = zeros(N,Global.M); y3 = zeros(N,Global.M);
    for i = 1: size(Off2,1)
        for j = 1 : Global.M
            [y1(i,j),~,MSE1(i,j)] = predictor(Off2(i,r1),Model1{j});
            [y2(i,j),~,MSE2(i,j)] = predictor(Off2(i,r2),Model2{j});
            [y3(i,j),~,MSE3(i,j)] = predictor(Off2(i,r3),Model3{j});
        end
    end
    y1 = [PopObj;y1]; y2=[PopObj;y2]; y3=[PopObj;y3];
    theta = Global.evaluated/Global.evaluation;
    PopObj = theta * y3 + (1 - theta) * (y1+y2)/2;  
    NV = size(V,1);
    yn = PopObj-repmat(min(PopObj,[],1),size(PopObj,1),1);
    Angle = acos(1-pdist2(yn,V,'cosine'));
    [~,associate] = min(Angle,[],2);
    Next = zeros(1,NV);
    for i = unique(associate)'
        current = find(associate==i);
        ED = sqrt(sum(yn(current,:).^2,2));
        [~,best] = min(ED);
        Next(i)  = current(best);
    end
    index = Next(Next~=0);
    Next = index;
    Parent = Off(Next,:);
    PopObj = PopObj(Next,:);
    V = V0.*repmat(max(PopObj,[],1)-min(PopObj,[],1),size(V0,1),1);
    Off2 = GA(Parent);
    Off = [Parent;Off2];
    g = g+1;
end
PopDec = Parent;
end

