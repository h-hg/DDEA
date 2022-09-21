%% Test functions' variable domain
function [Xmin, Xmax] = variable_domain(fname)
    switch fname
        case 'GRIEWANK'
            Xmin=-600; Xmax=600;
        case 'ACKLEY'
            Xmin=-32.768; Xmax=32.768;
        case 'ROSENBROCK'
            Xmin=-2.048; Xmax=2.048;
        case 'RASTRIGIN'
            Xmin=-5.12; Xmax=5.12;
        case 'ELLIPSOID'
            Xmin=-5.12; Xmax=5.12;
        case 'S_GRIEWANK'
            Xmin=-600; Xmax=600;
        case 'S_ACKLEY'
            Xmin=-32.768; Xmax=32.768;
        case 'S_ROSENBROCK'
            Xmin=-2.048; Xmax=2.048;
        case 'S_RASTRIGIN'
            Xmin=-5.12; Xmax=5.12;
        case 'S_ELLIPSOID'
            Xmin=-5.12; Xmax=5.12;
        case 'CEC05_F10'
            Xmin=-5; Xmax=5;
        case 'CEC05_F16'
            Xmin=-5; Xmax=5;
        case 'CEC05_F19'
            Xmin=-5; Xmax=5;
        otherwise
            Xmin=-5; Xmax=5;
    end
end
