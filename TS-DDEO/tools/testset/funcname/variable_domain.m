%------------------------------------------------------------------------
% This code is part of the program that produces the results in the following paper:
% Huixiang Zhen, Wenyin Gong, Ling Wang, Fei Ming, and Zuowen Liao. "Two-stage Data-driven Evolutionary Optimization for High-dimensional Expensive Problems", IEEE Transactions on Cybernetics, accepted, 2021.
% You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%----------------------------------------------------------------------------------------------------------------------------------------

function [Xmin, Xmax] = variable_domain(fname)
    switch fname
        case 'GRIEWANK'
            Xmin=-600; Xmax=600;
        case 'ACKLEY'
            Xmin=-32.768; Xmax=32.768;
        case 'ROSENBROCK'
            Xmin=-2.048; Xmax=2.048;
        case 'ELLIPSOID'
            Xmin=-5.12; Xmax=5.12;
        case 'CEC05_f10'
            Xmin=-5; Xmax=5;
        case 'CEC05_f16'
            Xmin=-5; Xmax=5;
        case 'CEC05_f19'
            Xmin=-5; Xmax=5;
        otherwise
            Xmin=-5; Xmax=5;
    end
end
