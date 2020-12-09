function [ POP,v ] = fly(POP,bu,bd,gbest,lbest,v,theta)
% Usage: [ POP,v ] = fly(POP,bu,bd,gbest,lbest,v,theta)
% -----------------------------------------------------------------
% Important Note: This code needs intralled SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:
% POP           - Population of Decision Variables
% gbest         - Global Best
% lbest         - Local Best
% bu            - Upper Boundary of c Decision Variables
% bd            - Lower Boundary of c Decision Variables
% v             - Velocity
% theta         - Generation Parameter
%
% Output: 
% POP           - Updated Population 
% v             - Updated Velocity
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
c1=1.49445;
c2=1.49445;
N=size(POP,1);
n=size(bu,2);
w=0.9-theta*0.5;
v=w*v+c1*rand(N,n).*(lbest(:,1:n)-POP(:,1:n))+c2*rand(N,n).*(ones(N,1)*gbest(1:n)-POP(:,1:n));
tvmax=0.5*ones(N,1)*(bu-bd);
I=find(v>tvmax);
v(I)=tvmax(I);
I=find(v<(0-tvmax));
v(I)=0-tvmax(I);
tPOP=POP(:,1:n)+v;

Tbd=ones(N,1)*bd;
Tbu=ones(N,1)*bu;
I=find(tPOP>Tbu);
tPOP(I)=Tbu(I)-(tPOP(I)-Tbu(I));
I=find(tPOP<Tbd);
tPOP(I)=Tbd(I)+(Tbd(I)-tPOP(I));

POP(:,1:n)=tPOP;
end