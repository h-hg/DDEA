function output = finess_ensemble_model_or(X)
global model_ensemble Dim

size_model_ensemble = length(model_ensemble);

set_f_RBF = [];
for i = 1:size_model_ensemble
    current_FUN = model_ensemble{i};
    f_RBF = current_FUN(X);
    set_f_RBF = [set_f_RBF, f_RBF];
end
f_mean = mean(set_f_RBF, 2);
for i = 1:size_model_ensemble
    diff(:, i) = (set_f_RBF(:, i) - f_mean).^2;
end
or = sum(diff, 2) / size_model_ensemble;

output = -or;