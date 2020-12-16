function [ Yhat,s ] = EnsemblePredicta( trees,S,nc,X )
N=size(nc,1);
Yhat=zeros(size(X,1),1);
s=zeros(size(X,1),1);
for i=1:size(X,1)
    f=zeros(N,1);
    for j=1:N
        f(j)=predict(trees{j,1},X(i,S(j,1:nc(j))));
    end
    Yhat(i)=mean(f);
    s(i)=std(f);
end

end