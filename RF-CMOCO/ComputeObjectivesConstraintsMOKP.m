function [ obj,Constraints ] =ComputeObjectivesConstraintsMOKP(POP,problem_name)
load(problem_name);
m=size(p,1);
c=size(p,2);
nc=size(w,1);
obj=[];
Constraints=[];
for i=1:m
    T=POP*p(i,:)';
    obj=[obj,T];
end
obj=0-obj;
for i=1:nc
    T=POP*w(i,:)'-cp(i);
    Constraints=[Constraints,T];
end
I=find(Constraints<0);
Constraints(I)=0;
end