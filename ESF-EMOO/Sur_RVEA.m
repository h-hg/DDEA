function [PopDec,PopObj] = Sur_RVEA(PopDec,PopObj,Ensemble,r,Global,V0)

Model1 = Ensemble{1};  Model2 = Ensemble{2};  Model3 = Ensemble{3};
r1 = r{1};  r2 = r{2};  r3 = r{3};
[f,~] = NDSort(PopObj,1);
fin = find(f==1);
rin = find(f~=1);
FObj = PopObj(fin,:);
FDec = PopDec(fin,:);
RDec = PopDec(rin,:);
RObj = PopObj(rin,:);

NDec = PopDec(fin,:);
NObj = FObj;
if length(fin) < 10
    NDec = [FDec;RDec(randperm(size(RDec,1),10-length(fin)),:)];
    NObj = [FObj;RObj(randperm(size(RDec,1),10-length(fin)),:)];
end
%  if length(fin) > 50
%     [ind,~]= kmeans(FObj,50);
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
%         NDec = FDec(Next,:);
%         NObj = FObj(Next,:);
% end
ind = randperm(size(NDec,1));
PopDec = NDec(ind,:);
PopObj = NObj(ind,:);

g = 1;  gmax = 20;  alpha = 2; V = V0;

while g <= gmax
    drawnow();
    OffDec = GA(PopDec);
    N = size(OffDec,1);
    y1 = zeros(N,Global.M); y2 = zeros(N,Global.M); y3 = zeros(N,Global.M); OffObj = zeros(N,Global.M); 
    for i = 1: N
        for j = 1 : Global.M
            [y1(i,j),~,MSE1(i,j)] = predictor(OffDec(i,r1),Model1{j});
            [y2(i,j),~,MSE2(i,j)] = predictor(OffDec(i,r2),Model2{j});
            [y3(i,j),~,MSE3(i,j)] = predictor(OffDec(i,r3),Model3{j});
        end
    end
    theta = Global.evaluated/Global.evaluation;
    OffObj = 0.5*(1-theta)*(y1+y2)/2 + 0.5*(1+theta)*y3;
%     OffObj =(y1+y2)/2;
    PopDec = [PopDec;OffDec];  PopObj = [PopObj;OffObj];
    index  = KEnvironmentalSelection(PopObj,V,(g/gmax)^alpha);
    PopDec = PopDec(index,:);
    PopObj = PopObj(index,:);
    % Adapt referece vectors
    if ~mod(g,ceil(gmax*0.1))
        V(1:Global.N,:) = V0.*repmat(max(PopObj,[],1)-min(PopObj,[],1),size(V0,1),1);
    end
    g = g + 1;
end

end

