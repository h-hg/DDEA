%%%*********************************************************************************************%%%
%% Implementation of Surrogate-assisted Hierarchical Particle Swarm Optimization (SHPSO)
%% H. Yu, Y. Tan, J. Zeng, C. Sun, Y. Jin, Surrogate-assisted hierarchical 
%% particle swarm optimization, Information Sciences, 454-455 (2018) 59-72.
%%%*********************************************************************************************%%%
%% This paper and this code should be referenced whenever they are used to 
%% generate results for the user's own research. 
%%%*********************************************************************************************%%%
%% This matlab code was written by Haibo Yu
%% Please refer with all questions, comments, bug reports, etc. to tyustyuhaibo@126.com
% 
%% Main process of SHPSO
% 
function [bestfit,hx,hf,fitcount,CE,gfs,pbest,pbestval,formula1,formula2]...
    = SHPSO_clean_version_for_share(fname,Dimension,Particle_Number,VRmin,VRmax,pos0,e0,hx,hf,fitcount,mf,CE,sn1,gfs)
% rand('seed',sum(100*clock));
ps=Particle_Number;                  %----种群规模
D=Dimension;                         %---维数
cc=[2.05 2.05];                      % acceleration constants  
iwt=0.7298;                          % 固定惯性权重
mv=0.5*(VRmax-VRmin);
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);            %---位置上界,防溢出
Vmin=repmat(-mv,ps,1);
Vmax=-Vmin;                          %---速度上界
vel=Vmin+2.*Vmax.*rand(ps,D);        % initialize the velocity of the particles
pos=pos0;e=e0;                       % initialize position and fitness
pbest=pos;  pbestval=e;              % initialize the pbest and the pbest's fitness value
[gbestval,gbestid]=min(pbestval);
gbest=pbest(gbestid,:);              % initialize the gbest and the gbest's fitness value
gbestrep=repmat(gbest,ps,1);

if D < 100
    gs=100;                          % < 100 dimension    建模训练样本数
elseif D == 100
    gs=150;                          % = 100 dimension
end

besty=1e200;
bestp=zeros(1,D);
formula1=0; formula2=0;              % 统计formula1和formula2的使用次数

gen = 0;

% Main loop
while fitcount < mf
    % Sample training samples
    [ghf,id]=sort(hf);               % 使用历史数据库中最优的若干个体构造全局代理模型
    ghf=ghf(1:gs);     ghx=hx(id(1:gs),:);

    %% R B F network
    ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
    % build global surrogate model
    spr=max(max(ghxd))/(D*gs)^(1/D);
    net=newrbe(ghx',ghf,spr);        % newrbe
    FUN=@(x) sim(net,x');           
    
    % Record the old best position and fitness
    besty_old=besty;
    bestp_old=bestp;
    % Find a near optimum of surrogate model by SLPSO
    maxgen=50*D; minerror=1e-6;
    [bestp,~] = SLPSO(D,maxgen,FUN,minerror,ghx);    % 利用SLPSO寻找全局模型的最优解
    % Exact evaluate the model optimum
    besty=feval(fname,bestp); 
    fitcount = fitcount + 1;
    if fitcount <= mf 
        CE(fitcount,:)=[fitcount,besty];
        if mod (fitcount,sn1)==0
            cs1=fitcount/sn1; gfs(1,cs1)=min(CE(1:fitcount,2));
        end
    end
% 
    % Record the new best position and fitness of surrogate
    besty_new=besty;
    bestp_new=bestp;
    % Update model optimum
    if besty_new < besty_old      
        besty=besty_new;    
        bestp=bestp_new;
        bestprep=repmat(bestp_new,ps,1);
    else
        besty=besty_old;
        bestp=bestp_old;
        bestprep=repmat(bestp_old,ps,1);
    end    
%         
    [~,ih,~]=intersect(hx,bestp,'rows');
    if isempty(ih)==1
        hx=[hx;bestp];  hf=[hf,besty];              % update history database 
    end
%     
    % Update the RBF model                                        
    if ghf(end) > besty                             % 选取最优样本集建立全局模型
        [ghf,id]=sort(hf);                          % 使用历史数据库中最优的若干个体构造全局代理模型
        ghf=ghf(1:gs); ghx=hx(id(1:gs),:);
        %% R B F network
%         ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));  
        spr=max(max(ghxd))/(D*gs)^(1/D);
        net=newrbe(ghx',ghf,spr);       % newrbe
        FUN=@(x) sim(net,x');   

    end    

    if besty < gbestval                             % 代理模型的近似最优小于当前代PSO种群的gbest,使用besty引导粒子向模型最优位置飞行
        [~,ip,~]=intersect(pbest,gbest,'rows');
        pbest(ip,:) =bestp;
        pbestval(ip)=besty;                         % 更新gbest对应的pbest
        gbestrep=bestprep;
        aa=cc(1).*rand(ps,D).*(pbest-pos)+cc(2).*rand(ps,D).*(gbestrep-pos);
        formula1=formula1+1;
    else
        aa=cc(1).*rand(ps,D).*(pbest-pos)+cc(2).*rand(ps,D).*(gbestrep-pos);
        formula2=formula2+1;
    end
    
    vel=iwt.*(vel+aa);                              % 带收缩因子的PSO算法
    vel=(vel>Vmax).*Vmax+(vel<=Vmax).*vel;
    vel=(vel<Vmin).*Vmin+(vel>=Vmin).*vel;
    pos=pos+vel;
    pos=((pos>=VRmin)&(pos<=VRmax)).*pos...
        +(pos<VRmin).*(VRmin+0.25.*(VRmax-VRmin).*rand(ps,D))+(pos>VRmax).*(VRmax-0.25.*(VRmax-VRmin).*rand(ps,D));
    
    %% Fitness estimation of new population
    e=FUN(pos);
    
    % Prescreening candidate solutions for exact evaluation
    candidx=find(e < pbestval);
    pos_trmem=pos(e < pbestval, :);                 % e-pbest strategy    
    [~,ih,ip]=intersect(hx,pos_trmem,'rows');
    if ~isempty(ip)==1        
        pos_trmem(ip,:)=[];
        e(candidx(ip))=hf(ih);
        candidx(ip)=[];
    end
    % Exact evaluate the prescreened candidates
    ssk=size(pos_trmem,1);
    e_trmem=zeros(1,ssk);        
    for k=1:ssk
        e_trmem(k)=feval(fname,pos_trmem(k,:));
        fitcount = fitcount + 1;
        if fitcount <= mf 
            CE(fitcount,:)=[fitcount,e_trmem(k)];
            if mod (fitcount,sn1)==0
                cs1=fitcount/sn1; gfs(1,cs1)=min(CE(1:fitcount,2));
            end
        end
        hx=[hx;pos_trmem(k,:)];hf=[hf,e_trmem(k)];        % update history database 
        kp=candidx(k);
        if e_trmem(k)<pbestval(kp)                        % update pbest
            pbest(kp,:)=pos_trmem(k,:);
            pbestval(kp)=e_trmem(k); 
        end
    end  
    
    % Update gbest
    [gbestval,tmp]=min(pbestval); 
    gbest=pbest(tmp,:);
    gbestrep=repmat(gbest,ps,1);                          % update the gbest    
    bestfit=min([gbestval,besty]);
    gen = gen + 1;
    
    fprintf(1,'Iteration: %d,  No.evaluation: %d,  Best: %e,  No.prescreen data: %d\n',gen,fitcount,bestfit,ssk);     
end
