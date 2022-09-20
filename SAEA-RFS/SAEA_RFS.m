addpath(genpath(pwd));
clear
warning off
global initial_flag
global lu
m = 10;  % population size
ts_num = 200;  % the number of initial training samples
K = 20;  % the number of sub-problems
sub_size = 100;  % maximum dimension of all sub-problems 
runnum = 25;  % the number of trial runs
nIter = 5; % maximum iterations of sub-problem optimization
results = zeros(15,runnum);
total_time = zeros(15,runnum);
for func_num = 1:1
    initial_flag = 0;
    %for overlapping function D = 905
    if func_num > 12 && func_num < 15
        D = 905; % dimensionality of the objective function.
    else
        D = 1000;
    end
    
    % Search Range
    if (func_num == 1 | func_num == 4 | func_num == 7 | func_num == 8 | func_num == 11 | func_num == 12 | func_num == 13 | func_num == 14 | func_num == 15 )
        lu = [-100*ones(1,D);100*ones(1,D)];
    end
    if (func_num == 2 | func_num == 5 | func_num == 9)
        lu = [-5*ones(1,D);5*ones(1,D)];
    end
    if (func_num == 3 | func_num == 6 | func_num == 10)
        lu = [-32*ones(1,D);32*ones(1,D)];
    end
    
    d = D;  % dimensions of the problem
    maxfe = 11*d;  % maxfe: maximal number of fitness evaluations
    
    % several runs
    for run = 1:runnum
        % generate initial archive
        XRRmin = repmat(lu(1, :), ts_num, 1);
        XRRmax = repmat(lu(2, :), ts_num ,1);
        archive = XRRmin + (XRRmax - XRRmin) .* lhsdesign(ts_num ,d);
        archive(:, d+1) = benchmark_func(archive(:,1:d), func_num);
        FES = ts_num;
        cFE = FES;
        cycle = 0;
        flag = 'cubic';
        tic;
        bestever = min(archive(:,d+1));
        fprintf('cycle = %d\tBest fitness: %e\n',cycle,bestever);
        array = bestever;
        %% Initial population
        XRRmin = repmat(lu(1, :), m, 1);
        XRRmax = repmat(lu(2, :), m ,1);
        Parent = XRRmin + (XRRmax - XRRmin) .* lhsdesign(m ,d);
        
        % main loop
        while (FES < maxfe)
            %% form sub-problem and training
            ts = archive;
            [lambda, gamma, dv , sp ] = RFS(sub_size,ts(:, 1:d),...
                ts(:, end), flag, K);
            %% Optimization
            [bestY, index] = min(archive(:, end));
            bestX = archive(index, 1:d);
            sub_eval_min = [];
            
            for i = 1:K
                sub_Parent = Parent(:,dv{i});   % the population of sub-problem
                fit_Parent = RBF_eval(sub_Parent,ts(sp{i}, dv{i}), ...
                        lambda{i}, gamma{i}, flag);
                    
                % use DE as a basic optimizer for each sub-problem
                for iter = 1:nIter
                    sub_child = DE(sub_Parent, bestX(1,dv{i}));
                    % Model evaluation
                    fit_child = RBF_eval(sub_child, ts(sp{i}, dv{i}), ...
                        lambda{i}, gamma{i}, flag);
                    fit_pop = [fit_Parent;fit_child];
                    [~,y] = sort(fit_pop);
                    com_pop = [sub_Parent;sub_child];
                    % Environment selection
                    sub_Parent = com_pop(y(1:m),:);
                    fit_Parent = fit_pop(y(1:m));
                end
             % save the best solution and its estimated value of sub-problem
                sub_eval_min(i,1) = fit_Parent(1);
                sub_eval_min(i,2:length(dv{i})+1) = sub_Parent(1,:); 
                
                % update the problem population
                Parent(:,dv{i}) = sub_Parent;  
            end
            
            %% Determine the solution to be evaluated using real funcion
            solution = bestX;
            [min_model,sign] =min(sub_eval_min(:,1));
            temp = RBF_eval(Parent(1,dv{sign}),ts(sp{sign},dv{sign}),lambda{sign},gamma{sign},flag);
            if temp <= min_model
                solution(:,dv{sign}) = Parent(1,dv{sign});
            else
                solution(:,dv{sign}) = sub_eval_min(sign,2:length(dv{sign})+1);
            end
            %% search difference
            diff = setdiff(solution,archive(:,1:d),'rows');
            diff(:,d+1) = benchmark_func(diff(:,1:d), func_num);
            FES = FES + 1;
            archive = [archive;diff];
            cFE = [cFE;FES];
            cycle = cycle+1;
            bestever = min(min(archive(:,d+1)),bestever);  % update the global solution
            array = [array;bestever];
            fprintf('cycle = %d\tBest fitness: %e\n',cycle,bestever);
        end
        list{func_num,run}(:,1) = cFE;
        list{func_num,run}(:,2) = array;
        fprintf('Run No.%d Done!\n', run);
        disp(['CPU time: ',num2str(toc)]);
        total_time(func_num,run) = toc;
        results(func_num, run) = bestever;
    end
end


