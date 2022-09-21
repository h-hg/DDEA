function output = finess_model_main(X)
global Dim model_RBF

f_RBF = model_RBF(X(:, 1:Dim));

output = f_RBF;