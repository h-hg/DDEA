function [ NPOP ] = SBX( POP,bu,bd,pc,n )
%SBX交叉
%入口参数：POP-父代，bu-上边界矩阵，bd-下边界矩阵，pc-交叉概率
%出口参数：NPOP-子代
NPOP=[];
eta_c=15;
N=size(POP,1);
C=size(bu,2);
y=1;
for i=1:n/2
    r1=rand;
    if r1<=pc
        A=randperm(N);
        k=i;
        y=A(1);
        if k==y 
            y=A(2);
        end
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

