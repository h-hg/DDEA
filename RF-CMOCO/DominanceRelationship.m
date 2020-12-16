function [x] = DominanceRelationship(a,b,m,c)
%Input:
%   a,b-individual
%   c-code length
%   m-objective number
%Output
%   x-dominance relationship
%       1-a dominates b
%       2-b dominates a
%       3-a = b
%       4-a and b are non-dominated by each other

t=0;
q=0;
p=0;
for i=1:m
    if a(1,c+i)<=b(1,c+i)
        t=t+1;
    end
    if a(1,c+i)>=b(1,c+i)
        q=q+1;
    end
    if a(1,c+i)==b(1,c+i)
        p=p+1;
    end
end
if t==m&p~=m
    x=1;
elseif q==m&p~=m
    x=2;
elseif p==m
    x=3;
else
    x=4;
end
end

