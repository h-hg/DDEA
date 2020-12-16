function [dis ] = MinDis(POP,train,c )
n=size(POP,1);
dis=zeros(n,1);
for i=1:n
    d=sum(abs(ones(size(train,1),1)*POP(i,1:c)-train(:,1:c)),2);
    dis(i)=min(d);
end
end

