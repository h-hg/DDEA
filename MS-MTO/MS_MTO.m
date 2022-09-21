function MS_MTO=MS_MTO(Tasks,pop,gen,rmp,reps,fit,achieve_num)

no_of_tasks=length(Tasks);
y_y=Tasks(no_of_tasks).dims+1;

minrange = Tasks(no_of_tasks).Lb(1:Tasks(no_of_tasks).dims);
maxrange = Tasks(no_of_tasks).Ub(1:Tasks(no_of_tasks).dims);
y=maxrange-minrange;

for rep=1:reps
    disp(rep);
    tic
    
    lhs=[];
    archive=[];
    FES=0;
    generation=0;
    mu=1;
    mum=1;
    distance_m=((fit-achieve_num)/2)/10;
    toall_distance=distance_m;
    
    %% Initial archive
    lhs = lhsdesign(achieve_num,Tasks(no_of_tasks).dims);
    for i=1:achieve_num
        lhs(i,1:Tasks(no_of_tasks).dims) = y.*lhs(i,1:Tasks(no_of_tasks).dims) + minrange; % decoding
        lhs(i,Tasks(no_of_tasks).dims+1)=Tasks(no_of_tasks).fnc(lhs(i,1:Tasks(no_of_tasks).dims));
    end
    archive=lhs;
    FES=achieve_num;
    generation=generation+1;
        
    [xxx,yy]=sort(archive(:,Tasks(no_of_tasks).dims+1));
    archive=archive(yy,:);
    best_one(rep,:)=archive(1,:);
    
    Toall_FES(rep,generation)=FES;
    Toall_BestFitss(rep,generation)=best_one(rep,Tasks(no_of_tasks).dims+1);
    
    disp(['MS-MTO Generation = ', num2str(generation), ' best factorial costs = ', num2str(best_one(rep,y_y))]);
    
    while FES<fit
       %% Initial or update RBF
        length_achieve=length(archive);
        flag1='cubic';
        %Global model
        [lambda1, gamma1]=RBF(archive(:,1:Tasks(no_of_tasks).dims),archive(:,Tasks(no_of_tasks).dims+1),flag1);
        Tasks(1).fnc=@(x) RBF_eval(x,archive(:,1:Tasks(no_of_tasks).dims),lambda1,gamma1,flag1);
        
        %Local model
        [lambda1, gamma1]=RBF(archive(1:achieve_num,1:Tasks(no_of_tasks).dims),archive(1:achieve_num,y_y),flag1);
        Tasks(2).fnc=@(x) RBF_eval(x,archive(1:achieve_num,1:Tasks(no_of_tasks).dims),lambda1,gamma1,flag1);
               
        %  select experienced individuals
        if length_achieve<pop/4
            length_pop=length_achieve;
            sub_pop=archive(1:length_pop,:);
        else
            length_pop=pop/4;
            num_best=1;
            num_choose=0.5*pop;
            if length_achieve<num_choose
                sub_pop(1:num_best,1:Tasks(no_of_tasks).dims)=archive(1:num_best,1:Tasks(no_of_tasks).dims);
                sub_pop1=archive(num_best+1:length_achieve,1:Tasks(no_of_tasks).dims);
                rndlist_=randperm(length_achieve-num_best);
                sub_pop(num_best+1:length_pop,1:Tasks(no_of_tasks).dims)=sub_pop1(rndlist_(1:length_pop-num_best),1:Tasks(no_of_tasks).dims);
            else
                sub_pop(1:num_best,1:Tasks(no_of_tasks).dims)=archive(1:num_best,1:Tasks(no_of_tasks).dims);
                sub_pop1=archive(num_best+1:num_choose,1:Tasks(no_of_tasks).dims);
                rndlist_=randperm(num_choose-num_best);
                sub_pop(num_best+1:length_pop,1:Tasks(no_of_tasks).dims)=sub_pop1(rndlist_(1:length_pop-num_best),1:Tasks(no_of_tasks).dims);
            end
        end
        % Normalization [0,1]
        sub_pop2(1:length_pop,1:Tasks(no_of_tasks).dims)= (sub_pop(1:length_pop,1:Tasks(no_of_tasks).dims) - minrange)./y;
        
        % The parameters linetype updating of SBX and Polynomial Mutation
        if generation>toall_distance
            toall_distance=toall_distance+distance_m;
            mu =mu+1;
            mum =mum+1;
        end
        
        % GMFEA optimization
        data_GMFEA=GMFEA(Tasks,pop,gen,rmp,sub_pop2,mu,mum);
        
        nvars=data_GMFEA.bestInd_data(1,1).rnvec;
        vars = y.*nvars + minrange; % decoding
        archive(length_achieve+1,1:Tasks(no_of_tasks).dims)=vars;
        archive(length_achieve+1,Tasks(no_of_tasks).dims+1)=Tasks(no_of_tasks).fnc(vars);
        
        nvars=data_GMFEA.bestInd_data(1,2).rnvec;
        vars = y.*nvars + minrange; % decoding
        archive(length_achieve+2,1:Tasks(no_of_tasks).dims)=vars;
        archive(length_achieve+2,Tasks(no_of_tasks).dims+1)=Tasks(no_of_tasks).fnc(vars);
        
        FES=FES+2;
        generation=generation+1;
        
       %% Update and sort the archive
        archive=unique(archive,'rows');
        [xxx,yy]=sort(archive(:,Tasks(no_of_tasks).dims+1));
        archive=archive(yy,:);
        best_one(rep,:)=archive(1,:);
        
        Toall_FES(rep,generation)=FES;
        Toall_BestFitss(rep,generation)=best_one(rep,Tasks(no_of_tasks).dims+1);
        
        disp(['MS-MTO Generation = ', num2str(generation), ' best factorial costs = ', num2str(best_one(rep,y_y))]);        
    end
    wall_clock_time(rep,1)=toc;
end
MS_MTO.wall_clock_time=wall_clock_time;
MS_MTO.best_one=best_one;
MS_MTO.Toall_FES=Toall_FES;
MS_MTO.Toall_BestFitss=Toall_BestFitss;
end