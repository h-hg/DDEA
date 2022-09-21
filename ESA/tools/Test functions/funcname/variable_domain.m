%%%*********************************************************************************************%%%
%% Test functions' variable domain
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

%% For 'FITNESS' 
% Xmin=-32.768;Xmax=32.768;%---Ackley 
% Xmin=-600;Xmax=600;%---Griewank
% Xmin=-2.048;Xmax=2.048;%---Rosenbrock
% Xmin=-5.12;Xmax=5.12;%---Ellipsoid 

%% For 'benchmark_func'
% Xmin=-5;Xmax=5; % cec05func 10
% Xmin=-5;Xmax=5; % cec05func 16
% Xmin=-5;Xmax=5; % cec05func 19
