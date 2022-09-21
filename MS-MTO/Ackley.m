function obj = Ackley(var)
%Ackley function
%   - var: design variable vector

    dim = length(var);
    sum1 = 0; sum2 = 0;
    for i = 1: dim
        sum1 = sum1 + var(i)*var(i);
        sum2 = sum2 + cos(2*pi*var(i));
    end
    avgsum1 = sum1/dim;
    avgsum2 = sum2/dim;    
    obj = -20*exp(-0.2*sqrt(avgsum1)) - exp(avgsum2) + 20 + exp(1);    
end

