function [ child ] = Mutation( p1,pm,Bd )
%mutation
%Input:
%   p1,p2-parents
%   Bd:constraint, 2-hospitals can be MTC,1- hospital can be TU.
%Output
%   child1,child2:child individual
c=size(p1,2);
child=p1;
if rand<=1
    for i=1:c
        if rand<pm
            switch child(i)
                case 0
                    if Bd(i)==2
                        if rand<0.5
                            child(i)=1;
                        else
                            child(i)=2;
                        end
                    else
                        child(i)=1;
                    end
                case 1
                    if Bd(i)==2
                        if rand<0.5
                            child(i)=0;
                        else
                            child(i)=2;
                        end
                    else
                        child(i)=0;
                    end
                case 2
                    if rand<0.5
                        child(i)=0;
                    else
                        child(i)=1;
                    end
            end
        end
    end
else
     k=randperm(c);
     i=2;
     while i<c & child(k(1))==child(k(i))
         i=i+1;
     end
     I=[k(1),k(i)];
     child(I(1:2))=child(I(2:-1:1));
     for i=1:2
        if Bd(I(i))<child(I(i))
            child(I(i))=Bd(I(i));
        end
     end
end


end

