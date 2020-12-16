function [NPOP]=mutation(POP,bu,bd,pm)
%变异
%入口参数：POP-种群，bu-变量上界，gmax-变量下界，pm-变异概率
%出口参数：NPOP变异后结果，
N=size(POP,1);
C=size(bu,2);
eta_m=15;
NPOP=POP(:,1:C);
for i=1:N
    k=i;
    NPOP(i,:)=POP(k,1:C);
    for j=1:C
        r1=rand;
        if r1<=pm
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

