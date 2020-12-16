function [ NPOP ] = MutationPopulation( POP,c,pm,Bd,N )
%mutation for the whole population
%Input:
%   POP-parent population
%   c-code length
%   Bd:constraint, 2-hospitals can be MTC,1- hospital can be TU.
%Output
%   NPOP-offspring population
n=size(POP,1);
NPOP=zeros(N,c);
for i=1:N
   NPOP(i,:) = Mutation( POP(unidrnd(n),1:c),pm,Bd );
end

end



