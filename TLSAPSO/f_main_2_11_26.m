function [evaluation_timestemp,evaluation_gbesttemp,max_itertemp,best_fit,hitcount_timestemp,global_approximate_output,local_approximate_output]=f_main()
clc
clear all
global dimension
global popsize
global max_length
global evaltimes %实际计算次数
global param1 %RBF1精度
global param2 %RBF2精度
global initial_flag
global global_approximate %全局估计次数
global local_approximate %局部估计次数

initial_flag=0;
dimension=30;
popsize=60;
%popsize=60;
% run_times=1;
max_length=10000;
c1=2.05;
c2=2.05;
w=0.7;
xmax=100;
xmin=-100;
vmax=100;
vmin=-100;
evaltimes=0; %评价次数
format short;
% 
% max_iteration=zeros(run_times,1);
% evaluation_times=[];
% evaluation_gbest=[];
% hitcount_times=[];

% for r=1:run_times
%     clc
    %clear all
    current_position=zeros(popsize,dimension);
    current_velocity=zeros(popsize,dimension);
    current_fitness=zeros(popsize,1);
    pbest=zeros(popsize,dimension);
    pbest_fitness=zeros(popsize,1);
    gbest=zeros(1,dimension);
    gbest_fitness=zeros(1,1);
    gbest_output=[];%zeros(max_length+1);
    evaltimes_output=[];%zeros(max_length+1); %输出评价次数
    evaltimes_output1=[];%zeros(max_length+1); %输出每一代的评价次数
    %hitCount=zeros(max_length,1); %记录碰撞次数
    whethercal=zeros(popsize,1);
    whethereval=zeros(popsize,1);
    archive=[]; %存放近若干代实际计算过的个体信息
    global_approximate_output=[];
    local_approximate_output=[];
    
    
    currpos_arch=zeros(popsize,dimension+3); %存储个体最近实际计算的位置信息
    nbthreshod=zeros(1,dimension); %每维上的阈值
    
    avggbest_fitness=0; %平均全局最优适应值
    % evaluation_times=[];
    % evaluation_gbest=[];
    % hitcount_times=[];
    global_flag=7;
    
    %  for r=1:run_times
    add_node=0;
    archive=[];
    sum_to_update=2000;
    currpos_arch=[];
    %初始化建立avi文件
    %      aviobj=avifile('research1.avi');
    hashtable=init_hash();%初始化一个Hash表
    hashtable_position=init_hash();
    hashtable_archive=init_hash();
    precise=4;
    totalHitCount=0;%到当前代为止个体与历史位置相碰撞的次数
    evaltimes=0;
    d=zeros(dimension,popsize+1);
    selected=zeros(dimension,popsize);
    %初始化过程
    for i=1:popsize
        for t=1:dimension
            current_velocity(i,t)=rand*(vmax-vmin)+vmin;
%             current_position(i,t)=rand*(xmax-xmin)+xmin;
        end
    end
    dxmin=xmin;
    dxmax=dxmin+(xmax-xmin)/popsize;
     
    for i=1:popsize
        for t=1:dimension
            %current_position(i,t)=rand*((dxmax-(dxmax-dxmin)/2)-(dxmin+(dxmax-dxmin)/2))+(dxmin+(dxmax-dxmin)/2);
            current_position(i,t)=rand*(dxmax-dxmin)+dxmin;
        end
        dxmin=dxmax;
        dxmax=dxmin+(xmax-xmin)/(popsize);
    end
%     plot(current_position(:,1),current_position(:,2),'o');

    
    for i=1:popsize
        %current_fitness(i,1)=fitness(current_position(i,:));
        temp_position=current_position(i,:);
        initial_flag=0;
        current_fitness(i,1)=benchmark_func(temp_position,global_flag);
        zput(hashtable_position,current_position(i,1:end),current_fitness(i,1),precise); %将实际计算的信息存入hash表中
        evaltimes=evaltimes+1;
        %存储初始化时个体实际计算的位置信息
        currpos_arch(i,1)=i;
        currpos_arch(i,2)=0;
        currpos_arch(i,3)=current_fitness(i,1);
        for t=1:dimension
            currpos_arch(i,t+3)=current_position(i,t);
        end
        archive(i,1)=i;%存储个体号
        archive(i,2)=0; %存储代数
        archive(i,3)=current_fitness(i,1); %存储适应值
        for t=1:dimension
            archive(i,t+3)=current_position(i,t);   %存储位置
        end
        zput(hashtable_archive,current_position(i,1:end),current_fitness(i,1),precise);
    end
    
    
    evaltimes_output(1)=evaltimes;
    evaltimes_output1(1)=evaltimes_output(1);
    global_approximate=0;
    local_approximate=0;
    global_approximate_output(1)=0;
        local_approximate_output(1)=0; 
        
    for i=1:popsize
        pbest(i,:)=current_position(i,:);
        pbest_fitness(i,:)=current_fitness(i,:);
    end
    
    gbest(1,:)=pbest(1,:);
    gbest_fitness(1,1)=pbest_fitness(1,1);
    for i=2:popsize
        if pbest_fitness(i,1)<gbest_fitness(1,1)
            gbest(1,:)=pbest(i,:);
            gbest_fitness(1,1)=pbest_fitness(i,1);
        end
    end
    gbest_output(1)=gbest_fitness(1,1); %输出初始化时最优值

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    totalupdate1=0; %统计net1样本更新个数
    %sum_to_update=1000;
    iter=1;
    fit_threshold=0.0001;%0.1时不稳定，
    fit_threshold1=0.0001;
    %while (evaltimes<=10000)&&(iter<=max_length)
    while (evaltimes<=10000) & (abs(gbest_fitness(1,1)+450)>=0.001)
        %fit_threshold=0.1-0.09*evaltimes/8000;
        currpos_archtemp=currpos_arch;
        for i=1:popsize
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---更新速度和位置
            for t=1:dimension
                %current_velocity(i,t)=w*current_velocity(i,t)+c1*rand*(pbest(i,t)-current_position(i,t))+c2*rand*(gbest(1,t)-current_position(i,t));
                current_velocity(i,t)=0.7298*(current_velocity(i,t)+c1*rand*(pbest(i,t)-current_position(i,t))+c2*rand*(gbest(1,t)-current_position(i,t)));
                if current_velocity(i,t)>vmax
                    current_velocity(i,t)=vmax;
                end
                if current_velocity(i,t)<vmin
                    current_velocity(i,t)=vmin;
                end
                current_position(i,t)=current_position(i,t)+current_velocity(i,t);
                if current_position(i,t)>xmax
                    current_position(i,t)=xmax;
                end
                if current_position(i,t)<xmin
                    current_position(i,t)=xmin;
                end
            end
        end
        xmax_t=zeros(dimension,1);
        xmin_t=zeros(dimension,1);
        for t=1:dimension
            xmax_t(t,1)=max(current_position(:,t));
            xmin_t(t,1)=min(current_position(:,t));
            nbthreshold(1,t)=abs(xmax_t(t,1)-xmin_t(t,1));
            xmax_t(t,1)=xmax_t(t,1)+0.25*(xmax_t(t,1)-xmin_t(t,1));%(xmax-xmin)/popsize;%
            xmax_t(t,1)=min(xmax,xmax_t(t,1));
            xmin_t(t,1)=xmin_t(t,1)-0.25*(xmax_t(t,1)-xmin_t(t,1));%(xmax-xmin)/popsize;%
            xmin_t(t,1)=max(xmin,xmin_t(t,1));
        end

        
        curr_index=1;
        for i1=1:size(currpos_arch,1)
            input_curr=1;
            for t=1:dimension
                if currpos_arch(i1,t+3)<xmin_t(t,1) || currpos_arch(i1,t+3)>xmax_t(t,1)
                    %                  if currpos_arch(i1,t+3)<t_min || currpos_arch(i1,t+3)>t_max
                    input_curr=0;
                    break;
                end
            end
            if input_curr==1
                for t=1:dimension
                    curr_pos(curr_index,t)=currpos_arch(i1,t+3);
                end
                curr_fitness(curr_index,1)=currpos_arch(i1,3);
                curr_index=curr_index+1;
            end
        end
        max_pos=zeros(dimension,1);
        min_pos=zeros(dimension,1);
        for t=1:dimension
            max_pos(t,1)=max(curr_pos(:,t));
            min_pos(t,1)=min(curr_pos(:,t));
        end
        %              maxmin_temp=max_pos-min_pos;
        maxmin_temp=abs(max_pos-min_pos);
        spread1=min(maxmin_temp(find(maxmin_temp>0)));
        %          spread1=min(max_pos-min_pos);
        sample_size=size(currpos_arch,1);
        add_node=20;
        timenet=newrb(curr_pos',curr_fitness',0.1,spread1,add_node,5);
        

        real_eval=0;
        for i=1:popsize
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---判断是否已经计算过，若未计算判断周围实际计算过的位置信息是否大于阈值，若是，构建局部模型对当前个体进行适应值估计，否
            %则用全局模型对当前个体进行适应值估计。
            val=zget(hashtable_position,current_position(i,1:end),precise); %从hash表中取值
            if isempty(val)
                current_fitness(i,1)=sim(timenet,current_position(i,:)');
                global_approximate=global_approximate+1;
                global_fitness(i,1)=current_fitness(i,1);
                whethereval(i,1)=0;
                whethercal(i,1)=0;
            else
                current_fitness(i,1)=val;%直接使用hash表中适应值，不需要计算
                currpos_arch(size(currpos_arch,1)+1,:)=[i iter current_fitness(i,1) current_position(i,:)];
                totalHitCount=totalHitCount+1;
                whethercal(i,1)=1;
                whethereval(i,1)=1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if whethereval(i,1)==0
                neiborornot=zeros(size(archive,1),1);
                for i1=1:size(archive,1)
                    neiborornot(i1,1)=0;
                    for t=1:dimension
                        %if abs(current_position(i,t)-archive(i1,3+t))>=nbthreshold_temp1
                            %if abs(current_position(i,t)-archive(i1,3+t))>=nbthreshold_temp2 
                            if (abs(current_position(i,t)-archive(i1,3+t))>=nbthreshold(1,t)/2) & (nbthreshold(1,t)>0)
                            neiborornot(i1,1)=1;
                            break;
                        end
                    end
                end
                nbsize=size(find(neiborornot==0),1);

                if nbsize>=dimension
                    neighbor_info=zeros(size(find(neiborornot==0),1),dimension);
                    neighbor_infofit=zeros(size(find(neiborornot==0),1),1);
                    indice=1;
                    for i1=1:size(archive,1)
                        if neiborornot(i1,1)==0
                            for t=1:dimension
                                neighbor_info(indice,t)=archive(i1,t+3);
                            end
                            neighbor_infofit(indice,1)=archive(i1,3);
                            indice=indice+1;
                        end
                    end
                    min_pos=zeros(dimension,1);
                    max_pos=zeros(dimension,1);
                    for t=1:dimension
                        min_pos(t,1)=min(neighbor_info(:,t));
                        max_pos(t,1)=max(neighbor_info(:,t));
                    end
                    maxmin_temp=abs(max_pos-min_pos);
                    max_min=min(maxmin_temp(find(maxmin_temp>0)));
                    spread2=max_min;
                    local_add_node=20;
                    spacenet=newrb(neighbor_info',neighbor_infofit',0.01,spread2,local_add_node,5); %局部代理模型
                    temp=current_fitness(i,1);
                    if(current_fitness(i,1)>sim(spacenet,current_position(i,:)'))
                        current_fitness(i,1)=sim(spacenet,current_position(i,:)');
                        local_approximate=local_approximate+1;
                    end
                    %                      current_fitness(i,1)=(current_fitness(i,1)+sim(spacenet,current_position(i,:)'))/2;
                    if current_fitness(i,1)<pbest_fitness(i,1)
                        real_eval=real_eval+1;
                        temp_position=current_position(i,:);
                        initial_flag=0;
                        current_fitness(i,1)=benchmark_func(temp_position,global_flag);
                        whethercal(i,1)=1;
                        difference=abs((temp-current_fitness(i,1))/current_fitness(i,1));
                        zput(hashtable_position,current_position(i,1:end),current_fitness(i,1),precise); %将实际计算的信息存入hash表中
                        if difference>fit_threshold || current_fitness(i,1)<pbest_fitness(i,1)%若差值大于给定阈值，则将实际计算个体的信息存入currpos_arch中
                            currpos_arch(size(currpos_arch,1)+1,:)=[i iter current_fitness(i,1) current_position(i,:)];
                            totalupdate1=totalupdate1+1;
                        end
                        evaltimes=evaltimes+1;
                        exist=0;
                        nbparam=1e-3;%0.0005+0.0005*evaltimes/8000;
                        for i1=1:size(archive,1)
                            disttoarchive=0;
                            for t=1:dimension
                                if abs(current_position(i,t)-archive(i1,t+3))>=nbparam%*(xmax-xmin)
                                    disttoarchive=1;
                                    break;
                                end
                            end
                            if disttoarchive==0
                                if abs(current_fitness(i,1)-archive(i1,3))/abs(current_fitness(i,1))>fit_threshold1
                                    disttoarchive=1;
                                end
                            end
                            if disttoarchive==0
                                exist=1;
                                break
                            end
                        end
                        if exist==0
                            recordnum=size(find(archive(:,1)==i),1);
                            if recordnum<10%30 %
                                archive_size=size(archive,1)+1;
                                archive(archive_size,1)=i;
                                archive(archive_size,2)=iter;
                                archive(archive_size,3)=current_fitness(i,1);
                                for t=1:dimension
                                    archive(archive_size,t+3)=current_position(i,t);
                                end
                            else
                                recordmatric=find(archive(:,1)==i);
                                particleminiter=min(archive(recordmatric,2));
                                replace=find(archive(:,1)==i & archive(:,2)==particleminiter);
                                archive(replace,1)=i;
                                archive(replace,2)=iter;
                                archive(replace,3)=current_fitness(i,1);
                                for t=1:dimension
                                    archive(replace,t+3)=current_position(i,t);
                                end
                            end
                        end
                        if current_fitness(i,1)<pbest_fitness(i,1)
                            pbest(i,:)=current_position(i,:);
                            pbest_fitness(i,1)=current_fitness(i,1);
                        end
                    end
                else
                    if current_fitness(i,1)<pbest_fitness(i,1)
                        real_eval=real_eval+1;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        temp=current_fitness(i,1);
                        %current_fitness(i,1)=fitness(current_position(i,:));
                        temp_position=current_position(i,:);
                        initial_flag=0;
                        current_fitness(i,1)=benchmark_func(temp_position,global_flag);
                        whethercal(i,1)=1;
                        difference=abs((temp-current_fitness(i,1))/current_fitness(i,1)); %计算差值
                        zput(hashtable_position,current_position(i,1:end),current_fitness(i,1),precise); %将实际计算的信息存入hash表中
                        if difference>fit_threshold  || current_fitness(i,1)<pbest_fitness(i,1)
                            currpos_arch(size(currpos_arch,1)+1,:)=[i iter current_fitness(i,1) current_position(i,:)];
                            totalupdate1=totalupdate1+1;
                        end
                        evaltimes=evaltimes+1;
                        exist=0;
                        nbparam=1e-3;%0.0005+0.0005*evaltimes/8000;
                        for i1=1:size(archive,1)
                            disttoarchive=0;
                            for t=1:dimension
                                if abs(current_position(i,t)-archive(i1,t+3))>=nbparam%*(xmax-xmin)
                                    disttoarchive=1;
                                    break;
                                end
                            end
                            if disttoarchive==0
                                if abs(current_fitness(i,1)-archive(i1,3))/abs(current_fitness(i,1))>fit_threshold1
                                    disttoarchive=1;
                                end
                            end
                            if disttoarchive==0
                                exist=1;
                                break
                            end
                        end
                        if exist==0
                            recordnum=size(find(archive(:,1)==i),1);
                            if recordnum<10 %
                                archive_size=size(archive,1)+1;
                                archive(archive_size,1)=i;
                                archive(archive_size,2)=iter;
                                archive(archive_size,3)=current_fitness(i,1);
                                for t=1:dimension
                                    archive(archive_size,t+3)=current_position(i,t);
                                end
                            else
                                recordmatric=find(archive(:,1)==i);
                                particleminiter=min(archive(recordmatric,2));
                                replace=find(archive(:,1)==i & archive(:,2)==particleminiter);
                                archive(replace,1)=i;
                                archive(replace,2)=iter;
                                archive(replace,3)=current_fitness(i,1);
                                for t=1:dimension
                                    archive(replace,t+3)=current_position(i,t);
                                end
                            end
                        end
                        if current_fitness(i,1)<pbest_fitness(i,1)
                            pbest(i,:)=current_position(i,:);
                            pbest_fitness(i,1)=current_fitness(i,1);
                        end
                    end
                end
            else
                exist=0;
                nbparam=1e-3%0.0005+0.0005*evaltimes/8000;
                for i1=1:size(archive,1)
                    disttoarchive=0;
                    for t=1:dimension
                        if abs(current_position(i,t)-archive(i1,t+3))>=nbparam%*(xmax-xmin)
                            disttoarchive=1;
                            break;
                        end
                    end
                    if disttoarchive==0
                        if abs(current_fitness(i,1)-archive(i1,3))/abs(current_fitness(i,1))>fit_threshold1
                            disttoarchive=1;
                        end
                    end
                    if disttoarchive==0
                        exist=1;
                        break
                    end
                end
                if exist==0
                    recordnum=size(find(archive(:,1)==i),1);
                    if recordnum<10 %
                        archive_size=size(archive,1)+1;
                        archive(archive_size,1)=i;
                        archive(archive_size,2)=iter;
                        archive(archive_size,3)=current_fitness(i,1);
                        for t=1:dimension
                            archive(archive_size,t+3)=current_position(i,t);
                        end
                    else
                        recordmatric=find(archive(:,1)==i);
                        particleminiter=min(archive(recordmatric,2));
                        replace=find(archive(:,1)==i & archive(:,2)==particleminiter);
                        archive(replace,1)=i;
                        archive(replace,2)=iter;
                        archive(replace,3)=current_fitness(i,1);
                        for t=1:dimension
                            archive(replace,t+3)=current_position(i,t);
                        end
                    end
                end
                if current_fitness(i,1)<pbest_fitness(i,1)
                    pbest(i,:)=current_position(i,:);
                    pbest_fitness(i,1)=current_fitness(i,1);
                end
                
            end
        end
        if real_eval==0
                for i=1:popsize
                temp=global_fitness(i,1);
                %current_fitness(individualflag(i),1)=fitness(current_position(individualflag(i),:));
                temp_position=current_position(i,:);
                initial_flag=0;
                current_fitness(i,1)=benchmark_func(temp_position,global_flag);
                %squaremse=squaremse+abs(log(temp-current_fitness(i,1)));
                difference=abs((temp-current_fitness(i,1))/current_fitness(i,1));
                zput(hashtable_position,current_position(i,1:end),current_fitness(i,1),precise); %将实际计算的信息存入hash表中
                if difference>fit_threshold || current_fitness(i,1)<pbest_fitness(i,1)
                    currpos_arch(size(currpos_arch,1)+1,:)=[i iter current_fitness(i,1) current_position(i,:)];
                    totalupdate1=totalupdate1+1;
                end
                evaltimes=evaltimes+1;
                exist=0;
                nbparam=1e-3;%0.0005+0.0005*evaltimes/8000;
                for i1=1:size(archive,1)
                    disttoarchive=0;
                    for t=1:dimension
                        if abs(current_position(i,t)-archive(i1,t+3))>=nbparam%*(xmax-xmin)
                            disttoarchive=1;
                            break;
                        end
                    end
                    if disttoarchive==0
                        if abs(current_fitness(i,1)-archive(i1,3))/abs(current_fitness(i,1))>fit_threshold1
                            disttoarchive=1;
                        end
                    end
                    if disttoarchive==0
                        exist=1;
                        break
                    end
                end
                if exist==0
                    recordnum=size(find(archive(:,1)==i),1);
                    if recordnum<10 %
                        archive_size=size(archive,1)+1;
                        archive(archive_size,1)=i;
                        archive(archive_size,2)=iter;
                        archive(archive_size,3)=current_fitness(i,1);
                        for t=1:dimension
                            archive(archive_size,t+3)=current_position(i,t);
                        end
                    else
                        recordmatric=find(archive(:,1)==i);
                        particleminiter=min(archive(recordmatric,2));
                        replace=find(archive(:,1)==i & archive(:,2)==particleminiter);
                        archive(replace,1)=i;
                        archive(replace,2)=iter;
                        archive(replace,3)=current_fitness(i,1);
                        for t=1:dimension
                            archive(replace,t+3)=current_position(i,t);
                        end
                    end
                    if current_fitness(i,1)<pbest_fitness(i,1)
                        pbest(i,:)=current_position(i,:);
                        pbest_fitness(i,1)=current_fitness(i,1);
                        %zput(hashtable,pbest(i,1:end),pbest_fitness(i,1),precise); %将实际计算的信息存入hash表中
                    end
                end
%                 end
            end
        end
        %end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---确定群体最优
        for i=1:popsize
            if pbest_fitness(i,1)<gbest_fitness(1,1)
                gbest(1,:)=pbest(i,:);
                gbest_fitness(1,1)=pbest_fitness(i,1);
            end
        end
        gbest_output(iter+1)=gbest_fitness(1,1);
        gbest_fitness(1,1)
        evaltimes_output(iter+1)=evaltimes;
        evaltimes_output1(iter+1)=evaltimes_output(iter+1)-evaltimes_output(iter);
        global_approximate_output(iter+1)=global_approximate;
        local_approximate_output(iter+1)=local_approximate;        
        iter=iter+1
        evaltimes
    end
    %iter-1
    %format long
    max_itertemp=iter-1
    for i=1:iter-1
        
        evaluation_timestemp(i)=evaltimes_output(i);
        evaluation_gbesttemp(i)=gbest_output(i);
        
    end
    best_fit=gbest_fitness(1,1);
    hitcount_timestemp=totalHitCount;

end