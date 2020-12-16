function [obj,Constraints ] = PredictObjectivesConstraintsMOKP(POP,Trees,S,NCtree,c,m,nc)

n=size(POP,1);
Constraints=zeros(n,nc);
obj=zeros(n,m);
for i=1:m
    [ obj(:,i),s ] = EnsemblePredicta( Trees(:,i),S(:,:,i),NCtree(:,i),POP(:,1:c) );
end

for i=m+1:m+nc
    [ Constraints(:,i-m),s ] = EnsemblePredicta( Trees(:,i),S(:,:,i),NCtree(:,i),POP(:,1:c) );
end
s=s;
end

