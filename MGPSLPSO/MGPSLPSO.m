function[gbestval,Itergbest,itereva,everyEVA,everyGBEST,time_output]  =  MGPSLPSO(Dimension,Max_Gen,VRmin,VRmax,problem)
 initial_time=cputime;
rng(sum(100*clock),'twister');
everyEVA=[];
everyGBEST=[];
maxEVA=Max_Gen;%z最大迭代次数

D=Dimension;

func_num=problem;

   


maxfe=maxEVA;

    M = 100;
    m = M + D/10;
%     c3 = D/M*0.01;
 c3 = 0;
    PL = zeros(m,1);

    for i = 1 : m
        PL(i) = (1 - (i - 1)/m)^log(sqrt(ceil(D/M)));%学习概率
    end


    %initialization
   ps=m;  

fiteva=zeros(ps,1);%每一代的对应的个体被实际计算时为1，否则为0 
 pastpos=[];
 pastposfit=[];
 time_output=[];
%%%%%%%%%
%%%%%%%%%拉丁超立方体初始化
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end

VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);

 v = zeros(m,D);

p=VRmin+(VRmax-VRmin).*lhsdesign(ps,D); %%%用拉丁超立方生成随机数
evacount=0;
    i2=1;


if D>30
    ID=3;%%%%%%%%%%ID是构建数据库的代数
else
    ID=2;
end

    %main loop
    while(evacount < maxfe)

      if i2<ID
         e(:,1)=evafit(p,func_num);
         fiteva(:,1)=1;
        fitest(:,1)=0;
        posfit=e;
        pastpos=[pastpos;p];
        pastposfit=[pastposfit;e];
        evacount=i2*ps;
        gbestval=min(pastposfit);
        Itergbest(i2,1)=gbestval;
        itereva(i2,1)=evacount; 
        everyEVA=[everyEVA;evacount];
        everyGBEST=[everyGBEST;gbestval];
        fprintf('Best fitness: %e\n', gbestval); 
      
        else
            GPfitness;
        Itergbest(i2,1)=gbestval;
        itereva(i2,1)=evacount; 
      end
        
        i2=i2+1;

        [posfit, gbestmid] = sort(posfit, 'descend');
        p = p(gbestmid,:);
        v = v(gbestmid,:);

        center = ones(m,1)*mean(p);

        randco1 = rand(m, D);

        randco2 = rand(m, D);

        randco3 = rand(m, D);
        winidxmask = repmat([1:m]', [1 D]);
        winidx = winidxmask + ceil(rand(m, D).*(m - winidxmask));

        pwin = p;
        for j = 1:D
                pwin(:,j) = p(winidx(:,j),j);
        end
        
         lpmask = repmat(rand(m,1) < PL, [1 D]);
         lpmask(m,:) = 0;
         v1 =  1*(randco1.*v + randco2.*(pwin - p) + c3*randco3.*(center - p));
         p1 =  p + v1;   
         
         
         v = lpmask.*v1 + (~lpmask).*v;         
         p = lpmask.*p1 + (~lpmask).*p;
         

            p = max(p, VRmin);
            p = min(p, VRmax);

    end

time_output = cputime - initial_time;
iternum=i2;
 fprintf('iternum: %d\n',iternum);

end