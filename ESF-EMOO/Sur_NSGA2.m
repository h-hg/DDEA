function [PopDec,PopObj] = Sur_NSGA2(PopDec,PopObj,Ensemble,r,Global)

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
% ind = randperm(size(PopDec,1));
% PopDec = PopDec(ind,:);
% PopObj = PopObj(ind,:);

NDec = PopDec(fin,:);
NObj = FObj;
if length(fin) < 10
    NDec = [FDec;RDec(randperm(size(RDec,1),10-length(fin)),:)];
    NObj = [FObj;RObj(randperm(size(RDec,1),10-length(fin)),:)];
end
ind = randperm(size(NDec,1));
PopDec = NDec(ind,:);
PopObj = NObj(ind,:);

g = 0;  gmax = 20;  N = size(PopDec,1);
[~,~,FrontNo,CrowdDis] = NSGA2_Selection(PopDec,PopObj,N);

while g <= gmax
    MatingPool = TournamentSelection(2,N,FrontNo,-CrowdDis);
    OffDec  = GA(PopDec(MatingPool,:));
    for i = 1: size(OffDec,1)
        for j = 1 : Global.M
            [y1(i,j),~,MSE1(i,j)] = predictor(OffDec(i,r1),Model1{j});
            [y2(i,j),~,MSE2(i,j)] = predictor(OffDec(i,r2),Model2{j});
            [y3(i,j),~,MSE3(i,j)] = predictor(OffDec(i,r3),Model3{j});
        end
    end
    theta = Global.evaluated/Global.evaluation;
    OffObj = 0.5*(1-theta)*(y1+y2)/2 + 0.5*(1+theta)*y3;
%     OffObj = theta * y3 + (1 - theta) * (y1+y2)/2;
%     theta = g/gmax;
%     OffObj = 0.5 * y3+ 0.5 * theta * y3 + 0.5 * (1 - theta) * (y1+y2)/2;
    PopDec = [PopDec;OffDec];  PopObj = [PopObj;OffObj];
    [PopDec,PopObj,FrontNo,CrowdDis] = NSGA2_Selection(PopDec,PopObj,N);
    g = g+1;
end
end

