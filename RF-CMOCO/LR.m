function [ Beta,Derta ] =LR( TrainData,Constraints,c,t )
nc=size(Constraints,2);
m=size(TrainData,2)-nc-c;
N=size(TrainData,1);
Beta=zeros(2,nc);
Derta=zeros(1,nc);

for i=1:nc
    x=Constraints(:,i);
    y=TrainData(:,c+m+i);
    I=find(y~=0);
    y(I)=1;
    B = glmfit(x, [y,ones(size(y,1),1)], 'binomial', 'link', 'logit');

    Beta(1,i)=B(1);
    Beta(2,i)=B(2);
    Derta(i)=(log(t/(1-t))-Beta(1,i))/Beta(2,i);
end

end

