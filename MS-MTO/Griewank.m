function obj = Griewank(var)
%GRIEWANK function
%   - var: design variable vector


    dim = length(var);
    sum1 = 0; sum2 = 1;
    for i = 1: dim
        sum1 = sum1 + var(i)*var(i);
        sum2 = sum2 * cos(var(i)/(sqrt(i)));
    end
    
    obj = 1+1/4000*sum1-sum2;
end