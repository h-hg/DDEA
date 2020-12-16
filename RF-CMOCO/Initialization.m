function [ pop ] = Initialization( n,c,Bd )
%Initialization of the poplation
%Input:
%   n-the number of Configurations
%   c-the number of hospitals
%   Bd:constraint, 2-hospitals can be MTC,1- hospital can be TU.
%Output
%   pop:the initial population

pop=zeros(n,c);
I1=find(Bd==1);
I2=find(Bd==2);
n1=n;
pop(1:n1,I1)=randint(n1,size(I1,2),[0,1]);
pop(1:n1,I2)=randint(n1,size(I2,2),[0,2]);


end

