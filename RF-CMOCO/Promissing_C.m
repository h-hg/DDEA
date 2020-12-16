function [ IR,ID ] = Promissing_C( SPOP,TPOP,c,m,nc)

n=size(SPOP,1);
np=zeros(n,1);
IR=[];
ID=[];
for i=1:size(SPOP,1)
    if sum(SPOP(i,c+m+1:c+m+nc))==0
        for j=1:size(TPOP,1)

            x=DominanceRelationship(SPOP(i,:),TPOP(j,:),m,c);
            if x==2|x==3
                np(i)=np(i)+1;
            end

        end
        if np(i)==0
            IR=[IR;i];
        else
            ID=[ID;i];
        end
    end
end
end

