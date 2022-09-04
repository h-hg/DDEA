%========================================================================================================================
%%%*********************************************************************************************************%%%
%% Implementation of Truncation-learning-driven surrogate assisted social learning particle swarm optimization for computationally expensive problem (TL-SSLPSO)
%% Haibo Yu , Li Kang , Ying Tan, Chaoli Sun , Jianchao Zeng,
%% Published in:Applied Soft Computing（Volume 97, Part A,?December 2020, 106812)
%%%*********************************************************************************************%%%
%% This paper and this code should be referenced whenever they are used to 
%% generate results for the user's own research. 
%%%*********************************************************************************************%%%
%% This matlab code was written by Haibo Yu
%% Please refer with all questions, comments, bug reports, etc. to tyustyuhaibo@126.com 
% % 
%=======================================================================================================================
clear,clc
rng('default');
rng('shuffle');
warning('off');
% d=51;% Stepped Cantilever Beam Design Problem
d=30; % d chosen from {30,50,100}
maxfes=1000;% temination by maximum number of FES

funcnum=1; % funcnum chosen from {1,2,,3,4,5,10,16,19}

if funcnum==1 %---Ellipsoid 
    fname='ellipsoid';
    Xmin=-5.12;Xmax=5.12;
    VRmin=-5.12;VRmax=5.12;
elseif funcnum==2 %---Rosenbrock
    fname='rosenbrock';
    Xmin=-2.048;Xmax=2.048;
    VRmin=-2.048;VRmax=2.048;
elseif funcnum==3 %---Ackley 
    fname='ackley';
    Xmin=-32.768;Xmax=32.768;
    VRmin=-32.768;VRmax=32.768;
elseif funcnum==4 %---Griewank
    fname='griewank';
    Xmin=-600;Xmax=600;
    VRmin=-600;VRmax=600;
elseif funcnum==5 %---Rastrigins 
    fname='rastrigin';
    Xmin=-5.12;Xmax=5.12;
    VRmin=-5.12;VRmax=5.12;
elseif funcnum==10 || funcnum==16 || funcnum==19 % CEC 2005 function F10/F16/F19
    fname='benchmark_func';
    Xmin=-5;Xmax=5;
    VRmin=-5;VRmax=5;
end

sn1=2;  gfs=zeros(1,fix(maxfes/sn1));   CE=zeros(maxfes,2);     % sampling setting according to FES

time_begin=tic;
runnum=30;

for run=1:runnum
    %---------------Initialization-----------------
    M = 100;      m = M + fix(d/10); % 种群规模 (下取整)
%%  Invariable C3
    c3 =0; % c3 >-1 guarantee convergence
% % %     beta=0.01; c3 = d/M*beta;     % epsilon    
%%
    PL = 1;
%     
    %initialization
    p = zeros(m, d); 
    v = zeros(m, d);    

    lu = [Xmin.* ones(1, d); Xmax.* ones(1, d)];
    FES = 0;    gen = 0;       
% % %     rand('seed', sum(100 * clock));

    % ----------*** Initialization *** Method -------- (IT1)
    XRRmin = repmat(lu(1, :), m, 1);
    XRRmax = repmat(lu(2, :), m, 1);
    p = XRRmin + (XRRmax - XRRmin) .* lhsdesign(m, d);          % 种群初始化: LHS
    fitness=zeros(1,m);
    for ii=1:m
        fitness(ii) = feval(fname,p(ii,:));
        FES=FES+1;
        if FES <= maxfes 
            CE(FES,:)=[FES,fitness(ii)];
            if mod (FES,sn1)==0
                cs1=FES/sn1;
                gfs(1,cs1)=min(CE(1:FES,2));
            end
        end
    end                                                     % 计算初始种群适应值
    hx=p;   hf=fitness;                                       % 初始化历史数据库
    
    [bestever,id] = min(fitness);
    gbpos=p(id,:);                                          % 精英保留策略: 种群最优个体位置
    %----------------End initialization------------------
%%
    % 通过设置参数group_num 来切换对比算法(RBF-SLPSO vs. TL-SSLPSO)
    group_num = 3; % // group_num=1 ---> RBF-SLPSO(S-SLPSO) // *** group_num>1--->Population division-based RBF-SLPSO (TL-SSLPSO) ***// take the value of 3 is of the best according to the sensitivity analysis
    %% 灵敏度参数分析：3，5，7，10，15
    %% 根据参数灵敏度实验结果，分成3组，整体效果较好
    
    %main loop
    while(FES < maxfes)                    
%         fprintf('Iteration: %d Fitness evaluation: %d Best fitness: %e\n', gen, FES, bestever);% 应该每评价一次，显示一次，而不是每迭代一次显示一次（应放在下面if后面）
        
        %% population segmentation
        if gen >= 1 && group_num~=1  % group_num=1 denotes the comparison algorithm: RBF-SLPSO algorithm
            % Equal division
            num_c=group_num;
            [~,idx]=sort(fitness); % 分割前，将数据按适应值从小到大进行排序
            idx=buffer(idx,ceil(m/num_c)); % idx 的每一列对应每一分割所包含样本的指标集      
            
            p_id=cell(num_c,1);% 各类中种群个体的指标集[初始化]
            pos_mean=cell(num_c,1);% 各类中种群个体集的平均位置
            centroid=cell(num_c,1);% 各类的质心位置
            nc=zeros(1,num_c);% 各类的类内样本数
            cf_mean=zeros(1,num_c);
            for i=1:num_c % divide database into m levels according to the fitness value
                % Equal division
                k=idx(:,i); k(k==0)=[];  
                pos_i=p(k,:);% positions of samples in ith cluster/group               
                cf=fitness(k);   cf_mean(i)=mean(cf); % mean fitness value of ith cluster/group
                centroid{i}=mean(pos_i); % the centroid of ith cluster/group
                [~,~,ip]=intersect(pos_i,p,'rows');
                if ~isempty(ip) == 1
                    p_id{i}=ip; % record the index of candidate individual belonging to the ith cluster/group
                    pos_mean{i}=mean(p(ip,:));
                else
                    p_id{i}=[];
                    pos_mean{i}=[];
                end              
                
                nc(i)=length(p_id{i});% number of samples in ith cluster/group
                
                % count the number of hx and p in ith cluster
                if i == 1 % the group of the highest-level sub-population
                    [~,ihh,~]=intersect(hx,pos_i,'rows');
                    if ~isempty(ihh)==1
                        num_hx=length(ihh);
                        num_p=nc(i)-num_hx;
                    else
                        num_hx=0;
                        num_p=nc(i);
                    end
                    num_app(gen,:)=[num_hx,num_p];
                end
                
            end

            [cf_mean,icf]=sort(cf_mean,'descend'); % 降序排列：according to fitness mean value of each cluster
            
            nc=nc(icf);                  % 排序后对应各类中样本数量
            p_id=p_id(icf);              % 排序后对应各类中个体的指标集
            pos_mean=pos_mean(icf);      % 排序后对应各类中种群个体的平均位置
            centroid=centroid(icf);      % 排序后对应各类样本的平均位置（质心）      
            
        end

       %% ********************* Phase 2 ************************ behavioral learning       
        if gen < 1
            %population sorting        
            [fitness,rank] = sort(fitness, 'descend');                 % 按降序（由大到小）进行排列   
            p = p(rank,:);
            v = v(rank,:); 
            %center position
            center = ones(m,1)*mean(p);
            %random matrix
            randco1 = rand(m, d);
            randco2 = rand(m, d);
            randco3 = rand(m, d);
            winidxmask = repmat([1:m]', [1 d]);
            winidx = winidxmask + ceil(rand(m, d).*(m - winidxmask));
%             winidx = m - floor(0.5*rand(m, d).*(m - winidxmask));    
            pwin = p;
            for j = 1:d
                pwin(:,j) = p(winidx(:,j),j);              % 每一个个体的第j维元素贡献学习向量的第j维
            end
            %learning
%             rand('seed', sum(100 * clock));
            lpmask = repmat(rand(m,1) < PL, [1 d]);
            lpmask(m,:) = 0;                              % 最优的第m个个体不学习
            v1 =  1*(randco1.*v + randco2.*(pwin - p) + c3*randco3.*(center - p));
            p1 =  p + v1;   
            v = lpmask.*v1 + (~lpmask).*v;                % 除最优个体外剩余个体都学习
            p = lpmask.*p1 + (~lpmask).*p;                % 保留最优个体到下一代：第m个个体不学习
        else
            for i = 1:group_num   
                p_index=p_id{i}(:,:);
                nn=length(p_index); % number of candidate inviduales in ith cluster/group                 
                if i == group_num % 排序后最优类种群个体的学习
                    if nn == 1 % 只包含一个种群个体
                        continue;% 最优类内的个体不学习
                    elseif group_num == 1 % *(comparison algorithm：RBF-SLPSO算法)* 自学习:实验发现，自学习容易破坏潜在最优解，误导种群的搜索
                        sam=[]; fit=[];                                
                        sam=p(p_index,:);
                        fit=fitness(p_index);
                        v_tmp=v(p_index,:);
                        [~,rank]=sort(fit,'descend');                    
                        sam=sam(rank,:);
                        fit=fit(rank);                                
                        %random matrix
                        randco1 = rand(nn, d);
                        randco2 = rand(nn, d);                    
                        randco3 = rand(nn, d);
                        winidxmask = repmat([1:nn]', [1 d]);
                        winidx = winidxmask + ceil(rand(nn, d).*(nn - winidxmask));   
                        pwin = sam;
                        for j = 1:d
                            pwin(:,j) = sam(winidx(:,j),j);              % 每一个个体的第j维元素贡献学习向量的第j维
                        end
                        %learning
                        lpmask = repmat(rand(nn,1) < PL, [1 d]);
                        lpmask(nn,:) = 0;                              % 最优的第m个个体不学习    
                        v1 =  1*(randco1.*v_tmp + randco2.*(pwin - sam) + c3*randco3.*(ones(nn,1)*pos_mean{i}(:,:) - sam)); % Main strategy:既向优秀个体学习保证算法的收敛性，同时又共享种群的平均信息保证算法的多样性,recall c3=0;
                        p1 =  sam + v1;   
                        v_tmp = lpmask.*v1 + (~lpmask).*v_tmp;                % 除最优个体外剩余个体都学习
                        sam = lpmask.*p1 + (~lpmask).*sam;                % 保留最优个体到下一代：第m个个体不学习

                        idd = p_index;
                        idd = idd(rank);
                        p(idd,:)=sam;
                        v(idd,:)=v_tmp;       

                        continue;  
                    else
                        continue
                    end
                    
                end                       
                    
                %% 需要找到比当前类层次高的所有类 (Demonstrator selection)
                id_demo=i+1:group_num; % 包含DB样本和种群个体的所有高层类 (Main strategy) 
                num_democluster=length(id_demo);                    

                for j = 1:d
                    %% competitive selection ---> RBF-SLPSO-division (TL-SSLPSO-B)
%                     if num_democluster == 1
%                         dim_cluster = num_c;
%                     else
%                         id=randperm(num_democluster);
%                         dim_cluster=id_demo(id(1:2)); % j 维度上 demonstrator *类指标*                        
%                         % 不同于level-based pso中比较高层类中随机对（random pair）个体的适应值,这里我们选取比较高层类间的适应值均值，选取较优类中的个体作为当前个体的demonstrators
%                         if cf_mean(dim_cluster(1)) < cf_mean(dim_cluster(2)) 
%                             dim_cluster=dim_cluster(1);
%                         else
%                             dim_cluster=dim_cluster(2);
%                         end                       
%                     end

                    %% random selection ---> RBF-SLPSO-division-random (TL-SSLPSO-R)
                    id=randi(num_democluster); % randomly select one---> for comparison
                    dim_cluster=id_demo(id);
                    %%                    
                    demo_p=p(p_id{dim_cluster}(:,:),:);
                    demonstrator=demo_p; % (Main strategy)：demonstrator太多，考虑选择性的保留若干代表性的历史样本
                    pwin(p_index,j) = demonstrator(randi(size(demonstrator,1),1,size(p_index,1)),j);
                end

                %learning
                lpmask = repmat(rand(nn,1) < PL, [1 d]); 
                %random matrix
                randco1 = rand(nn, d);
                randco2 = rand(nn, d);
                randco3 = rand(nn, d);     
                v1 =  1*(randco1.*v(p_index,:) + randco2.*(pwin(p_index,:)  - p(p_index,:)) + c3*randco3.*(ones(nn,1)*pos_mean{i}(:,:) - p(p_index,:))); % Main strategy

                p1 =  p(p_index,:) + v1;   
                v(p_index,:) = lpmask.*v1 + (~lpmask).* v(p_index,:);                % 除最优个体外剩余个体都学习
                p(p_index,:) = lpmask.*p1 + (~lpmask).* p(p_index,:);                % 保留最优个体到下一代：第m个个体不学习
            end
        end
        %boundary
        for i = 1:m
            p(i,:) = max(p(i,:), lu(1,:)); % 下界
            p(i,:) = min(p(i,:), lu(2,:)); % 上界
        end                       
        
        gen = gen + 1;  
        
% % %         % compute swarm diversity
% % %         p_mean=mean(p);
% % %         distance=real(sqrt(p_mean.^2*ones(size(p'))+ones(size(p_mean))*(p').^2-2*p_mean*(p')));
% % %         diversity(gen)=max(distance);
        
        %% RBF-modeling and evaluation
        % select training samples% 包围当前搜索种群的 DB 样本 ------> 局部模型（local）
        NS=2*(d+1);                                    % *** 算法性能参数 ***
        phdis=real(sqrt(p.^2*ones(size(hx'))+ones(size(p))*(hx').^2-2*p*(hx')));        
        [~,sidx]=sort(phdis,2);                        % 每行都进行排序   
        nidx=sidx; nidx(:,NS+1:end)=[];                % 个体的近邻样本指标集矩阵
        nid=unique(nidx);
        trainx=hx(nid,:);   
        trainf=hf(nid); 
        
        % radial basis function interpolation----(RBF-interpolation)
        flag='cubic';
        [lambda, gamma]=RBF(trainx,trainf',flag);
        FUN=@(x) RBF_eval(x,trainx,lambda,gamma,flag);

        fitness=FUN(p);     fitness=fitness';% model values
        
%% *** Expensive evaluation: individual-based evolution control **************** %%
        %% 计算适应度估值优于当前最优的个体的真实适应值，否则计算当前估值最优个体的真实适应值
        %% 实验发现，受限于代理模型精度及问题解空间特性的影响，算法的寻优过程（若干次迭代搜索）存在受困于某些真实的历史最优区域的风险。
        %% 尽管如此，算法结合分组学习机制，有能力跳出这些影响域,换言之，分组学习机制间接减小了算法陷入模型局部极值的概率，因为其保留了部分局部最优样本，这些样本不进行学习
        %% Two alternative infill sampling criteria
        pid = find(fitness < bestever);
        p_tmp = p(pid,:);
        f_tmp = fitness(pid);
        if ~isempty(p_tmp) == 1
            for i = 1 : size(p_tmp,1)
                [~,ih,~]=intersect(hx,p_tmp(i,:),'rows');
                if ~isempty(ih)==1
                    f_tmp(i)=hf(ih);
                else
                    f_tmp(i)=feval(fname,p_tmp(i,:));
                    FES=FES+1;
                    if FES <= maxfes
                        CE(FES,:)=[FES,f_tmp(i)];
                        if mod (FES,sn1)==0
                            cs1=FES/sn1;
                            gfs(1,cs1)=min(CE(1:FES,2));% the newly evaluated sample points are used to improve the accuracy of the surrogate model on one side, and on the other side, they are used as a promising leader in the evolution population
                        end
                    end
                    hx=[hx;p_tmp(i,:)];   hf=[hf,f_tmp(i)];                  % 更新历史数据库   
                end
                [bestever,ib] = min([f_tmp(i), bestever]);             % 更新全局最优         
                if ib==1
                    gbpos=p_tmp(i,:); 
                end
            end                
            fitness(pid) = f_tmp;
        else
            %% 计算适应度估值最优的个体的真实适应值
            [~,idx]=sort(fitness); 
            p_app=p(idx,:); f_app=fitness(idx);        
            [~,~,ip]=intersect(hx,p_app,'rows');
            p_app(ip,:)=[]; f_app(ip)=[]; % delete history samples
            
            sbest_pos=p_app(1,:);% performance ----> under limited computational cost, we prefer fast convergence performance                
            sbesty=feval(fname,sbest_pos);
            FES=FES+1;
            if FES <= maxfes
                CE(FES,:)=[FES,sbesty];
                if mod (FES,sn1)==0
                    cs1=FES/sn1;
                    gfs(1,cs1)=min(CE(1:FES,2));% the newly evaluated sample points are used to improve the accuracy of the surrogate model on one side, and on the other side, they are used as a promising leader in the evolutionary population
                end
            end
            hx=[hx;sbest_pos];   hf=[hf,sbesty];                  % 更新历史数据库   
            [bestever,ib] = min([sbesty, bestever]);             % 更新全局最优         
            if ib==1
                gbpos=sbest_pos; 
            end
            [~,ip,~]=intersect(p,p_app(1,:),'rows');
            fitness(ip)=sbesty;                
        end        
%% *****************************************************************         
 
    end
    
    fprintf('Runing Number: %d Fitness evaluation:%d Total generation: %d Best fitness:%e\n', run,FES,gen,bestever);
    
    gennum(run,:)=gen;
    
    gsamp1(run,:)=gfs;
end

best_samp=min(gsamp1(:,end));
worst_samp=max(gsamp1(:,end));
samp_mean=mean(gsamp1(:,end));
samp_median=median(gsamp1(:,end));
std_samp=std(gsamp1(:,end));
out1=[best_samp,worst_samp,samp_median,samp_mean,std_samp];
gsamp1_ave=mean(gsamp1,1);
gsamp1_log=log(gsamp1_ave);

for j=1:maxfes
    if mod(j,sn1)==0
        j1=j/sn1; gener_samp1(j1)=j;
    end
end

figure(1);
plot(gener_samp1,gsamp1_log,'.-r','Marker','.','Markersize',25,'LineWidth',2)
% semilogy(gener_samp1,gsamp1_ave,'.-r','Marker','.','Markersize',25,'LineWidth',2)
legend('TL-SSLPSO-R');
% xlim([100,maxfe]);
xlabel('Function Evaluation Calls ');
ylabel('Mean Fitness Value (Natural log)');
% % % title('2005 CEC Benchmark Function (F10)')
% % % title('2005 CEC Benchmark Function (F16)')
% % % title('2005 CEC Benchmark Function (F19)')
% % % title('Ackley Function')
% % % title('Griewank Function')
% % % title('Rastrigin Function')
% % % title('Rosenbrock Function')
% % % title('Ellipsoid Function')
set(gca,'FontSize',20);
time_cost=toc(time_begin);
time_ave=time_cost/runnum;