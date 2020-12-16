function [ NPOP ] = Cross( POP,TPOP,c,Bd )
%cross for the whole population

n=size(TPOP,1);
NPOP=zeros(2*n,c);
for i=1:n
    k=randint(1,1,[1,n]);
    [ NPOP(2*i-1,:),NPOP(2*i,:) ] = MultiCross(POP(i,1:c),TPOP(k,1:c),Bd);
end
end

