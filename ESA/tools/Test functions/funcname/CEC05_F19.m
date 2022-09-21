function f=CEC05_f19(x)
    global initial_flag
    initial_flag = 0;
    f=benchmark_func(x,19);
    f = f';
end