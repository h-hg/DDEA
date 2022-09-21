function f=CEC05_f16(x)
    global initial_flag
    initial_flag = 0;
    f=benchmark_func(x,16);
    f = f';
end