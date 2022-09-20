function V=mutation(X,bestX,F,mutationStrategy)
[NP,D]=size(X);

for i=1:NP
    %在[1 NP]中产生nrandI个互不相等的随机数，且与i皆不相等
    nrandI=5;
    r=randi([1,NP],1,nrandI);
    for j=1:nrandI
    equalr(j)=sum(r==r(j));   % 计数：产生的随机数中相同的个数
    end
    equali=sum(r==i);           % 行号与产生的随机数相同的次数
    equalval=sum(equalr)+equali;     % 相同的值
    while(equalval>nrandI)
        r=randi([1,NP],1,nrandI);
        for j=1:nrandI
        equalr(j)=sum(r==r(j));
        end
        equali=sum(r==i);
        equalval=sum(equalr)+equali;
    end
    
    switch mutationStrategy
        case 1
            %mutationStrategy=1:DE/rand/1;
            V(i,:)=X(r(1),:)+F*(X(r(2),:)-X(r(3),:));
        case 2
            %mutationStrategy=2:DE/best/1;
            V(i,:)=bestX+F*(X(r(1),:)-X(r(2),:));
        case 3
            %mutationStrategy=3:DE/rand-to-best/1;
            V(i,:)=X(i,:)+F*(bestX-X(i,:))+F*(X(r(1),:)-X(r(2),:));
        case 4
            %mutationStrategy=4:DE/best/2;
            V(i,:)=bestX+F*(X(r(1),:)-X(r(2),:))+F*(X(r(3),:)-X(r(4),:));
        case 5
            %mutationStrategy=5:DE/rand/2;
            V(i,:)=X(r(1),:)+F*(X(r(2),:)-X(r(3),:))+F*(X(r(4),:)-X(r(5),:));
        otherwise
            error('没有所指定的变异策略，请重新设定mutationStrategy的值');
    end
   
    
end
        
    