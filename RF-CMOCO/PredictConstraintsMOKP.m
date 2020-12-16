function [Constraints ] = PredictConstraintsMOKP(POP,Trees,S,NCtree,c,m,nc)

n=size(POP,1);
Constraints=zeros(n,nc);


for i=m+1:m+nc
    [ Constraints(:,i-m),s ] = EnsemblePredicta( Trees(:,i),S(:,:,i),NCtree(:,i),POP(:,1:c) );
end
end

