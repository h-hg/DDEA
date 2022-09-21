function f=CEC05_f10(x)
    global initial_flag
    initial_flag = 0;
    f=benchmark_func(x,10);
    f = f';
end