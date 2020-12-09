function [ obj ] =compute_objectives(POP,c,problem_name)
% Usage: [ obj ] =compute_objectives(POP,c,problem_name)
% -----------------------------------------------------------------
% Important Note: This code needs intralled SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:
% problem_name  - Benchmark Problem
% c             -No. of Decision Variables
% POP           -Population of Decision Variables
%
% Output: 
% obj           - Calculated Objective Value
%
    %%%%    Authors:    Handing Wang, Yaochu Jin, John Doherty
    %%%%    University of Surrey, UK
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage
    %%%%    DATE:       May 2018
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:
%Handing Wang, Yaochu Jin, John Doherty, Committee-based Active Learning for Surrogate-Assisted Particle Swarm Optimization of Expensive Problems, IEEE Transactions on Cybernetics, vol.47, no.9, pp.2664-2677, 2017.
%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
obj=[];
n=size(POP,1);
switch problem_name
    case 'Ackley'
        a=20;
        b=0.2;
        cc=2*pi;
        obj=0-a*exp(0-b*(mean(POP(:,1:c).^2,2)).^0.5)-exp(mean(cos(cc*POP(:,1:c)),2))+a+exp(1);
    case 'Rastrigin'
        obj=10*c+sum(POP(:,1:c).^2-10*cos(2*pi*POP(:,1:c)),2);
    case 'Rosenbrock'
        P=100*(POP(:,2:c)-POP(:,1:c-1).^2).^2;
        obj=sum(P,2)+sum((POP(:,1:c-1)-1).^2,2);
    case 'Griewank'
        P=[1:c].^0.5;
        P=ones(size(POP,1),1)*P;
        obj=sum(POP(:,1:c).^2,2)/4000+1-prod(cos(POP(:,1:c)./P),2);
    case 'Ellipsoid'
        P=[1:c];
        P=ones(size(POP,1),1)*P;
        obj=sum(P.*(POP(:,1:c).^2),2);     
end


end

