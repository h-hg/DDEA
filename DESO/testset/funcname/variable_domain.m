function [Xmin, Xmax] = variable_domain(fname)
    switch fname
        case 'GRIEWANK'
            Xmin=-600; Xmax=600;
        case 'ACKLEY'
            Xmin=-32.768; Xmax=32.768;
        case 'ROSENBROCK'
            Xmin=-2.048; Xmax=2.048;
        case 'Ellipsoid'
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
