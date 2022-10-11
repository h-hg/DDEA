function [YDec,obj] = Infill_Solutions_Selection(Dec,Obj,PopDec,PopObj,V0)

[Dec1,in] = unique(Dec,'rows');
Obj = Obj(in,:);
lia = ismember(Dec1,PopDec,'rows');
if sum(~lia)~=0
    Dec = Dec1(~lia,:);
    Obj = Obj(~lia,:);
end

[f,~] = NDSort([PopObj;Obj],1);
fin = f == 1;
n = size(Obj,1);
if sum(fin(end-n+1:end)) ~=0
    Dec = Dec(fin(end-n+1:end),:);
    Obj = Obj(fin(end-n+1:end),:);
end

V = V0.*repmat(max(Obj,[],1)-min(Obj,[],1),size(V0,1),1);
Obj = Obj-repmat(min(Obj,[],1),size(Obj,1),1);
l = randperm(size(V,1),5);
lamda  = V(l,:);
dis = pdist2(Obj,lamda,'cosine');
[mind,id] = min(dis,[],1);
YDec = Dec(id,:);
obj = Obj(id,:);
end

