function U=crossover(X,V,CR,crossStrategy)
[NP,Dim]=size(X);
global lu
switch crossStrategy
    %crossStrategy=1:binomial crossover
    case 1
        for i=1:NP
            jRand=floor(rand*Dim);%由于jRand要在[0,1)*Dim中取值，故而用floor
            for j=1:Dim
                if rand<CR||j==jRand
                    U(i,j)=V(i,j);
                else
                    U(i,j)=X(i,j);
                end     
            end    
        end
    %crossStrategy=2:Exponential crossover
    case 2
        for i=1:NP
            j=floor(rand*Dim);%由于j在[0,1)*Dim中取值，故而用floor
            L=0;
            U=X;
            U(i,j)=V(i,j);
            j=mod(j+1,D);
            L=L+1;
            while(rand<CR&&L<Dim)
                U(i,j)=V(i,j);
                j=mod(j+1,D);
                L=L+1;
            end
        end
        
    otherwise
        error('没有所指定的交叉策略，请重新设定crossStrategy的值'); 
end
%% 越界重置
for i=1:NP
    for j=1:Dim
        while (U(i,j)>lu(2,j)) || (U(i,j)<lu(1, j))
            U(i,j)=lu(1,j)+rand*(lu(2,j) - lu(1,j));
        end
    end

end
end
        
