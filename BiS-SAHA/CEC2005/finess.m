 function Y = finess(X)

global func_num
f = benchmark_func(X,func_num);
Y = [X, f];


