function [ NPOP ] = Nondominated_C( POP,c,m,nc)
%the non-dominated solution with constraints
%Input:
%   m-objective number
%   c-code length
%   POP-population
%   nc-the number of constrains
%Output
%   NPOP-selected individuals
n=size(POP,1);
Constraint =POP(:,c+m+1:c+m+nc);
nc=size(Constraint,2);

for i=1:nc
    if max(Constraint(:,i))~=0
        Constraint(:,i)=Constraint(:,i)/max(Constraint(:,i));
    end
end
Con=sum(Constraint,2);

np=zeros(n,1);
NPOP=[];
for i=1:n
    for j=1:n
        if j~=i
            x=DominanceRelationship_C(POP(i,:),POP(j,:),Con(i),Con(j),m,c);
            if x==2
                np(i)=np(i)+1;
            end
        end
    end
    if np(i)==0
        NPOP=[NPOP;POP(i,:)];
    end
end
end

