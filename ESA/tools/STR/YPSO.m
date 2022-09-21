function [xbest,fbest] = YPSO(name,Range)
c1=2.05;c2=2.05;
phi = c1 + c2;
format long;
global Dim
%------初始化种群的个体------------
N=50;
M=100;
Seed=lhsdesign(N,Dim);
v=lhsdesign(N,Dim);
N=length(Seed(:,1));
D=length(Seed(1,:));

for i=1:D
    Seed(:,i)=Seed(:,i)*(Range(i,2)-Range(i,1))+Range(i,1);
    v(:,i)=v(:,i)*(Range(i,2)-Range(i,1))+Range(i,1);
end

%------先计算各个粒子的适应度，并初始化Pi和Pg----------------------


p(:,1)=feval(name,Seed);

y=Seed;


pg = Seed(N,:);             %Pg为全局最优
fpg=p(N,1);

for i=1:(N-1)

    if p(i,1)<fpg

        pg=Seed(i,:);

    end

end

%------进入主要循环，按照公式依次迭代------------

for t=1:M

    for i=1:N
        ksi = 2*0.729 / abs(2 - phi - sqrt(phi^2 - 4*phi));
        v(i,:) = v(i,:)+c1*rand*(y(i,:)-Seed(i,:))+c2*rand*(pg-Seed(i,:));
        v(i,:) = ksi*v(i,:);
        Seed(i,:)=Seed(i,:)+v(i,:);
%%%%%%%%%%%%%%%%%%%判断位置变动是否在变量的设计域范围内
        for j=1:D
            if Seed(i,j)<Range(j,1)
                Seed(i,j)=Range(j,1);
            end
            if Seed(i,j)>Range(j,2)
                Seed(i,j)=Range(j,2);
            end
            if v(i,j)<Range(j,1)
                v(i,j)=Range(j,1);
            end
            if v(i,j)>Range(j,2)
                v(i,j)=Range(j,2);
            end
        end
    end
    fseed=feval(name,Seed);  
    for i=1:N
        
        if fseed(i)<p(i,1)

            p(i,1)=fseed(i);

            y(i,:)=Seed(i,:);

        end

        if p(i,1)<fpg

            pg=y(i,:);
            fpg=p(i,1);

        end

    end
end
xbest=pg;
fbest=fpg;
end



