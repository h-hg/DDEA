function output = finess_ensemble_model(X)
global model_ensemble Dim

X = X(:, 1:Dim);
size_model_ensemble = length(model_ensemble);
set_f_RBF = [];
for i = 1:size_model_ensemble
    current_FUN = model_ensemble{i};
    f_RBF = current_FUN(X);
    set_f_RBF = [set_f_RBF, f_RBF];
end

f_ensemble = max(set_f_RBF, [], 2);


output = f_ensemble;