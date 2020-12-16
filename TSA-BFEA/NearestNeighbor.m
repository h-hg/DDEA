function [I ] = NearestNeighbor(POP,train,c )
n=size(POP,1);
I=zeros(n,1);
for i=1:n
    d=sum(abs(ones(size(train,1),1)*POP(i,1:c)-train(:,1:c)),2);
    [A,I(i)]=min(d);
end
I=unique(I);
end

