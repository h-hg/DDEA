function [objective,rnvec,funcCount] = fnceval(Task,rnvec)
%FNCEVAL function: evaluate function in the unified search space
    d = Task.dims;
    nvars = rnvec(1:d);
    minrange = Task.Lb(1:d);
    maxrange = Task.Ub(1:d);
    y=maxrange-minrange;
    vars = y.*nvars + minrange; % decoding
    
    objective=Task.fnc(vars);
    funcCount=1;
end