function [ POP ] = DeleteSame( POP,c )
%Delete the same individuals in NSGA-II
%Input:
%   POP-parent population
%   c-code length
%Output
%   NPOP-population after sorted
i=1;
while i<=size(POP,1)-1
    I=[];
    for j=i+1:size(POP,1)
        Is=find(POP(i,1:c)==POP(j,1:c));
        if size(Is,2)==c
            I=[I;j];
        end
    end
    POP(I,:)=[];
    i=i+1;
end


end

