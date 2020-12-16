function [ Bd ] = Constraint (c)
%Generate the constraint of hispitals
%Input:
%   c-the number of hospitals
%Output
%   Bd:constraint, 2-hospitals can be MTC,1- hospital can be TU.
Bd=ones(1,c);
Bd(1)=2;
Bd(2)=2;
Bd(4)=2;
Bd(7)=2;

end

