function out = function_run(in,func_id,global_flag)
    if func_id <= 8
        out = fitness(in,func_id);
    else if (func_id > 8 && func_id <= 12)
            out = benchmark_func(in,global_flag);
        else if (func_id == 21 || func_id <= 22)
                out = benchmark_func(in,global_flag);
            else
                out = fitness(in, func_id);
            end
        end
    end
end