function [NPOP]=mutation(POP,bu,bd,pm,n)
% Usage: [NPOP]=mutation(POP,bu,bd,pm,n)
%
% Input:
% bu            -Upper Bound
% bd            -Lower Bound
% POP           -Input Population
% pm            -Mutation Probability
% n             -Population Size
%
% Output: 
% NPOP          -Output Population
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
N=size(POP,1);
C=size(bu,2);
eta_m=15;
NPOP=POP(:,1:C);
for i=1:n
%     k=randperm(N);
%     k=k(1);
    k=i;
    NPOP(i,:)=POP(k,1:C);
    for j=1:C
        r1=rand;
        if r1<=pm %±äÒì¸ÅÂÊ
            y=POP(k,j);
            yd=bd(j);yu=bu(j);
            if y>yd
                if (y-yd)<(yu-y)
                    delta=(y-yd)/(yu-yd);
                else
                    delta=(yu-y)/(yu-yd);
                end
                r2=rand;
                indi=1/(eta_m+1);
                if r2<=0.5
                    xy=1-delta;
                    val=2*r2+(1-2*r2)*(xy^(eta_m+1));
                    deltaq=val^indi-1;
                else
                    xy=1-delta;
                    val=2*(1-r2)+2*(r2-0.5)*(xy^(eta_m+1));
                    deltaq=1-val^indi;
                end
                y=y+deltaq*(yu-yd);
                NPOP(i,j)=min(y,yu);NPOP(i,j)=max(y,yd);
            else
                NPOP(i,j)=rand*(yu-yd)+yd;
            end
        end
    end
end
end

