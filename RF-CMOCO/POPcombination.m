function [ POP ] = POPcombination( POP1,POP2,c )
POP=POP1;
for i=1:size(POP2,1)
    T=sum(abs(POP(:,1:c)-ones(size(POP,1),1)*POP2(i,1:c)),2);
    if isempty(find(T==0))
        POP=[POP;POP2(i,:)];
    end
end


end

