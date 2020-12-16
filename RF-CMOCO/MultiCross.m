function [ child1,child2 ] = MultiCross(p1,p2,Bd)
%Multi point cross
%Input:
%   p1,p2-parents
%   Bd:constraint, 2-hospitals can be MTC,1- hospital can be TU.
%Output
%   child1,child2:child individual

Id=find(p1~=p2);
if size(Id,2)~=0
    child1=Mutation( p1,0.5,Bd );
    child2=Mutation( p1,0.5,Bd );    
elseif size(Id,2)>=0.3*size(p1,2)
    child1=p1;
    child2=p2;
    k=randint(1,1,[1,size(Id,2)]);
    Ii=randperm(size(Id,2));
    Ii=Ii(1:k);
    I=Id(Ii);
    child1(I)=p2(I);
    child2(I)=p1(I);
else
    Id=find(p1==p2);
    child1=p1;
    child2=p2;
    for i=1:size(Id,2)
        switch child1(Id(i))
            case 0
                if Bd(Id(i))==2
                    if rand<0.5
                        child1(Id(i))=1;
                        child2(Id(i))=2;
                    else
                        child1(Id(i))=2;
                        child2(Id(i))=1;
                    end
                else
                    child1(Id(i))=1;
                end
            case 1
                if Bd(Id(i))==2
                    if rand<0.5
                        child1(Id(i))=0;
                        child2(Id(i))=2;
                    else
                        child1(Id(i))=2;
                        child2(Id(i))=0;
                    end
                else
                     child1(Id(i))=0;
                end
            case 2
                if rand<0.5
                    child1(Id(i))=0;
                    child2(Id(i))=1;
                else
                    child1(Id(i))=1;
                    child2(Id(i))=0;
                end
        end
    end
end
end

