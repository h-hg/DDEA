function [ Trees,S,NCtree ] = SurrogateObjConRF(TrainData,c,m,nc,ntree )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Trees=[];
S=zeros(ntree,c,m+nc);
NCtree=zeros(ntree,m+nc);
for i=1:m
     [trees,S(:,:,i),NCtree(:,i)]=GenerateRFa( TrainData(:,1:c),TrainData(:,c+i),ntree);
     Trees=[Trees,trees];
end

for i=m+1:m+nc
    [trees,S(:,:,i),NCtree(:,i)] = GenerateRFa( TrainData(:,1:c),TrainData(:,c+i),ntree);
    Trees=[Trees,trees];
end

end

