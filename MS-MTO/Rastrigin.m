function obj = Rastrigin(var)
%Rastrigin function
%   - var: design variable vector

    dim = length(var);
    obj = 10*dim;
    for i=1:dim
        obj=obj+(var(i)^2 - 10*(cos(2*pi*var(i))));
    end
end