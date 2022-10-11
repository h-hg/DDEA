function [ NPOP ] = SBX( POP,bu,bd,pc,n )
% Usage: [ NPOP ] = SBX( POP,bu,bd,pc,n )
%
% Input:
% bu            -Upper Bound
% bd            -Lower Bound
% POP           -Input Population
% pc            -Crossover Probability
% n             -Population Scale
%
% Output: 
% NPOP          -Output Population with 2n Solutions
%
    %%%%    Authors:    Handing Wang, Yaochu Jin, Chaoli Sun, John Doherty
    %%%%    University of Surrey, UK and Taiyuan University of Science and Technology, China.
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage
    %%%%    DATE:       May 2018
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:

%Handing Wang, Yaochu Jin, Chaoli Sun, John Doherty, Offline data-driven evolutionary optimization using selective surrogate ensembles, IEEE Transactions on Evolutionary Computation, Accepted.

%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
NPOP=[];
eta_c=15;
N=size(POP,1); %种群大小
C=size(bu,2); %决策变量维数
y=1;
for i=1:n
    r1=rand; %随机概率
    if r1<=pc %交叉概率
        A=randperm(N); %打乱
        k=i;

        if A(2)<A(1)
            y=A(2);
        else
            y=A(1);
        end
        if k==y
            k=A(3); %k与y保证不一样
        end
        d=(sum((POP(y,1:C)-POP(k,1:C)).^2)).^0.5; %欧式距离
        if k~=y
            for j=1:C
                par1=POP(y,j);par2=POP(k,j);
                yd=bd(j);yu=bu(j);
                r2=rand;
                if r2<=0.5
                    y1=min(par1,par2);y2=max(par1,par2);
                    if (y1-yd)>(yu-y2)
                        beta=1+2*(yu-y2)/(y2-y1);
                    else
                        beta=1+2*(y1-yd)/(y2-y1);
                    end
                    expp=eta_c+1;beta=1/beta;alpha=2.0-beta^(expp);
                    r3=rand;
                    if r3<=1/alpha
                        alpha=alpha*r3;expp=1/(eta_c+1.0);
                        betaq=alpha^(expp);
                    else
                        alpha=1/(2.0-alpha*r3);expp=1/(eta_c+1);
                        betaq=alpha^(expp);
                    end
                    chld1=0.5*((y1+y2)-betaq*(y2-y1));
                    chld2=0.5*((y1+y2)+betaq*(y2-y1));   
                    aa=max(chld1,yd);
                    bb=max(chld2,yd);
                    if rand>0.5
                        NPOP(2*i-1,j)=min(aa,yu);
                        NPOP(2*i,j)=min(bb,yu);
                    else
                        NPOP(2*i,j)=min(aa,yu);
                        NPOP(2*i-1,j)=min(bb,yu);
                    end
                else
                    NPOP(2*i-1,j)=par1;
                    NPOP(2*i,j)=par2;
                end
            end
        end
    end
    
end
end

