function [ obj ] =compute_objective(POP,c,problem_name)
% Usage: [ obj ] =compute_objectives(POP,c,problem_name)

% Input:
% problem_name  - Benchmark Problem
% c             -No. of Decision Variables
% POP           -Population of Decision Variables
%
% Output:
% obj           - Calculated Objective Value
%

obj=[];
n=size(POP,1);
switch problem_name
    case 'Ellipsoid'
        P=[1:c];
        P=ones(size(POP,1),1)*P;
        obj=sum(P.*(POP(:,1:c).^2),2);
    case 'Rosenbrock'
        P=100*(POP(:,2:c)-POP(:,1:c-1).^2).^2;
        obj=sum(P,2)+sum((POP(:,1:c-1)-1).^2,2);
    case 'Rastrigin'
        obj=10*c+sum(POP(:,1:c).^2-10*cos(2*pi*POP(:,1:c)),2);
    case 'Ackley'
        a=20;
        b=0.2;
        cc=2*pi;
        obj=0-a*exp(0-b*(mean(POP(:,1:c).^2,2)).^0.5)-exp(mean(cos(cc*POP(:,1:c)),2))+a+exp(1);
    case 'Griewank'
        P=[1:c].^0.5;
        P=ones(size(POP,1),1)*P;
        obj=sum(POP(:,1:c).^2,2)/4000+1-prod(cos(POP(:,1:c)./P),2);
    case 'rastrigin_rot_func'
        obj = fun(POP(:,1:c),10);
    case 'hybrid_rot_func1'
        obj = fun(POP(:,1:c),16);
    case 'hybrid_rot_func2_narrow'
        obj = fun(POP(:,1:c),19);
end
end

