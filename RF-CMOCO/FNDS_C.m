function [ NPOP ] = FNDS_C( POP,c,m,nc)
%fast non-dominated sort by Deb with constraints handling
%Input:
%   POP-parent population
%   c-code length
%   m-objective number
%   nc-constraints number
%Output
%   NPOP-population after sorted
n=size(POP,1);
Constraint =POP(:,c+m+1:c+m+nc);
%normalization for constrains
for i=1:nc
    if max(Constraint(:,i))~=0
        Constraint(:,i)=Constraint(:,i)/max(Constraint(:,i));
    end
end
Con=sum(Constraint,2);
S=zeros(n,n);
np=zeros(n,1);
POP=[POP,zeros(n,1)];
NPOP=[];
F=[];
rank=1;
for i=1:n
    sp=[];
    for j=1:n
        if j~=i
            x=DominanceRelationship_C(POP(i,:),POP(j,:),Con(i),Con(j),m,c);
            if x==1
               sp=[sp,j];
            elseif x==2
                np(i)=np(i)+1;
            end
        end
    end
    if np(i)==0
        POP(i,m+c+nc+1)=1;
        F=[F,i];
    end
    S(i,1:size(sp,2))=sp;
end
while size(F,2)~=0
    Q=[];
    for i=1:size(F,2)
        p=F(i);
        j=1;
        while j<=n
            if S(p,j)~=0;
                k=S(p,j);
                np(k)=np(k)-1;
                if np(k)==0
                    POP(k,m+c+nc+1)=rank+1;
                    Q=[Q,k];
                end
            end
            j=j+1;
        end
    end
    rank=rank+1;
    F=Q;
end
NPOP=POP;
end


