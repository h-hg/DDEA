function output = BiS_SAHA
% [1] Zhihai Ren, Chaoli Sun, Ying Tan, Guochen Zhang, Shufen Qin. A 
% bi-stage surrogate assisted hybrid algorithm for expensive optimization 
% problems[J]. Complex & Intelligent Systems, 2021.
cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));

global initial_flag
global Dim Size_pop D
global func_num
global DataBase original_DB
global lb ub a b
global best_solution
global num_eval
global minerror
global model_RBF model_ensemble
global flag Nd

initial_flag = 0;
[a, b] = read_boundary;

lb = a * ones(1, Dim);
ub = b * ones(1, Dim);
MaxEval = 11 * Dim;
if Dim > 30
    MaxEval = 1000;
end

Size_DataBase = round(MaxEval*5/11);
Size_pop = Size_DataBase;
DataBase = a + (b-a) * lhsdesign(Size_DataBase, Dim);
DataBase = finess(DataBase);
original_DB = DataBase;
DataBase = sortrows(DataBase, Dim + 1);
best_solution = DataBase(1,:);
num_eval = 0;
num_eval = num_eval + Size_pop;

m = round(MaxEval*5/11); d = Dim;
pop_t = DataBase(1:m, :);
PL = ones(m,1);
c3 = 0;
Xmin = lb(1, 1);
Xmax = ub(1, 1);
lu = [Xmin* ones(1, d); Xmax* ones(1, d)];
fitness = pop_t(:, Dim + 1)';
p = pop_t(1:m, 1:Dim);
hx = p;
hf = fitness;
[~,idx]=sort(fitness); 
p_app=p(idx,:); f_app=fitness(idx);
[~,~,ip]=intersect(hx,p_app,'rows');
p_app(ip,:)=[];
f_app(ip)=[];

[~,idx]=sort(hf);   idx=idx(1:m);
p=hx(idx,:);        fitness=hf(idx);  
v = zeros(m, d); 
p;    fitness;  

eval_1 = [best_solution, num_eval];
eval_2 = [best_solution, num_eval];
eval_3 = [best_solution, num_eval];
eval_4 = [best_solution, num_eval];
eval_5 = [best_solution, num_eval];
D = Dim;
maxgen=50*Dim; minerror=1e-6;

F=0.8;%scaling factor
CR=0.8;%crossover rate
flag = 0;
flag_pop = 0;
AA = 0;

gen = 1;       
while 1

    if num_eval >= MaxEval
        break;
    end

    sub_DataBase = DataBase;  
    model_ensemble = build_ensemble_model(sub_DataBase);
    fitness = finess_ensemble_model(p);
    %% SLPSO
    num_gen = 100;
    while mod(gen, num_gen) ~= 0           
        [fitness,rank] = sort(fitness, 'descend');
        p = p(rank,:);
        v = v(rank,:); 
        center = ones(m,1)*mean(p);
        randco1 = rand(m, d);
        randco2 = rand(m, d);
        randco3 = rand(m, d);
        winidxmask = repmat([1:m]', [1 d]);
        winidx = winidxmask + ceil(rand(m, d).*(m - winidxmask));
        pwin = p;
        for j = 1:d
            pwin(:,j) = p(winidx(:,j),j);
        end
        lpmask = repmat(rand(m,1) < PL, [1 d]);
        lpmask(m,:) = 0;
        v1 =  1*(randco1.*v + randco2.*(pwin - p) + c3*randco3.*(center - p));
        p1 =  p + v1;   
        v = lpmask.*v1 + (~lpmask).*v;
        p = lpmask.*p1 + (~lpmask).*p;
        for i = 1:m
            p(i,:) = max(p(i,:), lu(1,:));
            p(i,:) = min(p(i,:), lu(2,:));
        end
        fitness = finess_ensemble_model(p);
        fitness=fitness';
        U = [p, fitness'];
        gen = gen + 1;
    end
    
    if mod(gen, num_gen) == 0
        gen = gen + 1;
        flag_pop = 1;
        flag = 1;
    end
    if flag_pop
        flag_pop = 0;
        sub_DataBase = init_pop_by_Clustering(DataBase);
        fitness = sub_DataBase(:, Dim + 1)';
        p = sub_DataBase(:, 1:Dim);
        fprintf('gen = %d eval = %d best = %e\n', gen, num_eval, best_solution(end, Dim + 1));
    end

    if 1
        sub_DataBase = DataBase;
        model_ensemble = build_ensemble_model(sub_DataBase);
        [A, B] = min(U(:, Dim + 1));
        eval_individual = U(B, :);
        eval_5 = [eval_5; [eval_individual, num_eval]];
        [~,~,ip]=intersect(DataBase(:, 1:Dim), U(:, 1:Dim),'rows');
        U(ip,:)=[];
        if ~isempty(U)
            if eval_individual(1, Dim + 1) < best_solution(1, Dim + 1)
                predictor_value = finess_ensemble_model_or(U(:, 1:Dim));
                U = [U(:, 1:Dim), predictor_value];
                U = sortrows(U, Dim + 1);
                eval_individual = U(1, :);
                eval_individual = finess(eval_individual(1, 1:Dim));
                Update_DB(eval_individual);
                eval_1 = [eval_1; [eval_individual, num_eval]];
                if num_eval >= MaxEval
                    break;
                end
            end
        end
        AA = AA + 1;
        fprintf('gen = %d eval = %d best = %e\n', gen, num_eval, best_solution(end, Dim + 1));
    end

    if num_eval >= round(MaxEval*6/11);
        sub_DataBase = DataBase;  
        [NP_sub_DataBase, ~] = size(sub_DataBase);
        randIndex = randperm(NP_sub_DataBase);
        S = sub_DataBase(randIndex, 1:Dim); Y = sub_DataBase(randIndex, Dim + 1); flag = 'cubic';
        [lambda, gamma]=RBF(S,Y,flag);
        model_RBF = @(x) RBF_eval(x,S,lambda,gamma,flag);
        gs_2 = round(MaxEval*0.5/11);
        gs = gs_2;
        sub_DataBase = sort_dis(DataBase);
        sub_DataBase = sub_DataBase(1:gs, :); 
        [sub_lb, sub_ub] = get_boundary(sub_DataBase);
        U = Multiple_iterations(sub_lb, sub_ub);
        eval_individual = U;
        if eval_individual(1, Dim + 1) < best_solution(1, Dim + 1)
            eval_individual = finess(eval_individual(1, 1:Dim));
            Update_DB(eval_individual);
            eval_3 = [eval_3; [eval_individual, num_eval]];
        end
        fprintf('gen = %d eval = %d best = %e\n', gen, num_eval, best_solution(end, Dim + 1));
    end
    flag = 0;
end
output = best_solution(end, Dim + 1);
filename = [pwd, '/', num2str(func_num), '_', num2str(Dim), '_', num2str(Nd), '.mat'];
save(filename, 'best_solution');

