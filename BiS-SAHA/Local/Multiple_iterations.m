function U = Multiple_iterations_main(sub_lb, sub_ub)

global Dim Size_pop
global minerror

pop = lhsdesign(Size_pop, Dim);
for i = 1:Dim
    pop(:, i) = sub_lb(1, i) + (sub_ub(1, i) - sub_lb(1, i)) * pop(:, i);
end
predictor_value = finess_model_main(pop);
pop = [pop(:, 1:Dim), predictor_value];
[~, num_best] = min(pop(:, Dim + 1));
bestX = pop(num_best, 1:(Dim + 1));
pop_t = pop;

F=0.8;%scaling factor
CR=0.8;%crossover rate
maxIteration = 150;
Generation = 0;
flag_b = 0;
while Generation  <= maxIteration
    V = mutation(pop_t(:, 1:Dim), F, bestX(end, 1:Dim));
    U = crossover(pop_t(:, 1:Dim), V, CR);
    U = Limited_boundary(U, sub_lb, sub_ub);
    predictor_value = finess_model_main(U);
    U = [U(:, 1:Dim), predictor_value];
    for i= 1:Size_pop
        if U(i, Dim + 1) < pop_t(i, Dim + 1)
            pop_t(i, 1:Dim + 1) = U(i, 1:Dim + 1);
        end
    end
    [~, Num_Xg] = min(pop_t(:, Dim + 1));
    bestX = [bestX; pop_t(Num_Xg, 1:Dim + 1)];
    Generation = Generation + 1;

    if abs(bestX(end, Dim + 1) - bestX(end - 1, Dim + 1)) < minerror
        flag_b = flag_b + 1;
    else
        flag_b = 0;
    end
    if flag_b == 20
        break;
    end
end

% fit = finess(bestX(:, 1:Dim));
U = pop_t(Num_Xg, :);

