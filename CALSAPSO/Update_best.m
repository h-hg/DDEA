function [ POP,gbest,lbest] = Update_best(POP,gbest,lbest)
% Usage: [ POP,gbest,lbest] = Update_best(POP,gbest,lbest)
% -----------------------------------------------------------------
% Important Note: This code needs intralled SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:
% POP           - Population of Decision Variables
% gbest         - Global Best
% lbest         - Local Best
%
% Output: 
% POP           - Updated Population 
% gbest         - Updated Global Best
% lbest         - Updated Local Best
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
[best,Ib]=min(POP(:,end));
if best<=gbest(end)
    gbest=POP(Ib,:);
end
I=find(POP(:,end)<=lbest(:,end));
lbest(I,:)=POP(I,:);
end