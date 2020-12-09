function [NPOP]= selectData( POP,n,c,f)
% Usage: [NPOP]= selectData( POP,n,c,f)
% -----------------------------------------------------------------
% Important Note: This code needs intralled SURROGATE TOOLBOX(https://sites.google.com/site/srgtstoolbox/)
% -----------------------------------------------------------------
% Input:
% n             -No. of Selected Data Points
% c             -No. of Decision Variables
% POP           -Population of Decision Variables
% f             -Norm
%
% Output: 
% NPOP          - Selected Data
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
NPOP=[];
I=1;
NPOP=[NPOP;POP(I,:)];
POP(I,:)=[];
while size(NPOP,1)<n & size(POP,1)>0
    d=zeros(size(POP,1),1);
    for i=1:size(POP,1)
        dis=sum((abs(NPOP(:,1:c)-ones(size(NPOP,1),1)*POP(i,1:c))).^f,2).^(1/f);
        d(i)=min(dis);
    end
    [A,I]=max(d);
    NPOP=[NPOP;POP(I,:)];
end
end