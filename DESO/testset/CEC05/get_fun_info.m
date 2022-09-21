function [L,U,opt_f,err] = get_fun_info(fun,D) 
% Get the lower-, upper-bound and the optimal value for the function 'fun' 
 
% Get the global optimal value of fun 
switch(fun) % fun from 1 to 25 are CEC'05 functions 
    case {1,2,3,4}, 
        LB = -100; UB = 100; 
        opt_f = -450; 
    case {5}, 
        LB = -100; UB = 100; 
        opt_f = -310;  
    case {6}, 
        LB = -100; UB = 100; 
        opt_f = 390; 
    case {7}, 
        LB = 0; UB = 600; 
        opt_f = -180; 
    case {8}, 
        LB = -32; UB = 32; 
        opt_f = -140; 
    case {9}, 
        LB = -5; UB = 5; 
        opt_f = -330; 
    case {10}, 
        LB = -5; UB = 5; 
        opt_f = -330; 
    case {11}, 
        LB = -0.5; UB = 0.5; 
        opt_f = 90; 
    case {12}, 
        LB = -100; UB = 100; 
        opt_f = -460; 
    case {13}, 
        LB = -3; UB = 1; 
        opt_f = -130; 
    case {14}, 
        LB = -100; UB = 100; 
        opt_f = -300; 
    case {15}, 
        LB = -5; UB = 5; 
        opt_f = 120; 
     case {16}, 
        LB = -5; UB = 5; 
        opt_f = 120; 
     case {17}, 
        LB = -5; UB = 5; 
        opt_f = 120; 
     case {18}, 
        LB = -5; UB = 5; 
        opt_f = 10; 
     case {19}, 
        LB = -5; UB = 5; 
        opt_f = 10; 
     case {20}, 
        LB = -5; UB = 5; 
        opt_f = 10;  
     case {21}, 
        LB = -5; UB = 5; 
        opt_f = 360; 
     case {22}, 
        LB = -5; UB = 5; 
        opt_f = 360; 
     case {23}, 
        LB = -5; UB = 5; 
        opt_f = 360; 
     case {24}, 
        LB = -5; UB = 5; 
        opt_f = 260; 
     case {25}, 
        LB = -2; UB = 5; 
        opt_f = 260; 
    case {26}, % Sphere function 
        LB = -100; UB = 100; 
        opt_f = 0; 
    case {27}, % Rastrigin 
        LB = -5.12; UB = 5.12; 
        opt_f = 0; 
    case {28}, %Six Hump Camel bsck 
        LB = -5; UB = 5; 
        opt_f = -1.031628; 
    case {29}, %Step 
        LB = -100; UB = 100; 
        opt_f = 0; 
    case {30}, %Rosenbrock 
        LB = -2; UB = 2; 
        opt_f = 0; 
    case {31}, %Ackley 
        LB = -32; UB = 32; 
       opt_f = 0; 
    case {32}, %Griewank 
        LB = -600; UB = 600; 
        opt_f = 0; 
    case {33}, %Salomon 
        LB = -100; UB = 100; 
        opt_f = 0; 
    case {34}, %Normalized Schwefel 
        LB = -512; UB = 512; 
        opt_f = -418.9828872724338; 
    case {35}, %Quartic 
        LB = -1.28; UB = 1.28; 
        opt_f = 0; 
    case {36}, %Rotated hyper-ellipsoid 
        LB = -100; UB = 100; 
        opt_f = 0; 
    case {37}, %Norwegian function 
        LB = -1.1; UB = 1.1; 
        opt_f = 1; 
    case {38}, %Alpine 
        LB = -10; UB = 10; 
        opt_f = 0; 
    case {39}, %Branin 
        LB = -5; UB = 15; 
        opt_f = 0.397887; 
    case {40}, %Easom 
        LB = -100; UB = 100; 
        opt_f = -1; 
    case {41}, % Goldstein and Price 
        LB = -2; UB = 2; 
        opt_f = 3; 
    case {42}, % Shubert 
        LB = -10; UB = 10; 
        opt_f = -186.7309; 
    case {43}, % Hartmann 
        LB = 0; UB = 1; 
        opt_f = -3.86278; 
    case {44}, % Shekel 
        LB = 0; UB = 10; 
        opt_f = -10.5364; 
    case {45}, % Levy 
        LB = -10; UB = 10; 
        opt_f = 0; 
    case {46}, % Michalewicz 
        LB = 0; UB = pi; 
        opt_f = -9.66015; 
    case {47}, % Shifted Griwank 
        LB = -600; UB = 600; 
        opt_f = -180; 
    case {48}, % Design of a Gear Train 
        LB = 12; 
    	UB = 60; 
        opt_f = 2.7e-12; 
    case {49}, % Pressure Vessel 
        LB = [1.125 0.625 0 0]; 
        UB = [12.5 12.5 240 240]; 
        opt_f = 7197.72893; 
    case {50}, % Tripod 
        LB = -100; %[-100 -100]; 
        UB = 100; % [100 100]; 
        opt_f = 0; 
    case {51}, % Compression Spring 
        LB = [1 0.6 0.207]; 
        UB = [70 3 0.5]; 
        opt_f = 2.6254214578; 
end 
 
% If LB and UB are not vectors make them vectors 
sl = size(LB); 
 
if (sl(1)*sl(2) == 1) % LB and UB are scalers 
    L = LB*ones(1,D); 
    U = UB*ones(1,D); 
else 
    L = LB; 
    U = UB; 
end 
 
% Admissible error 
 
if fun < 6 
    err = 1e-6; 
elseif fun >= 6 && fun < 17 
    err = 1e-2; 
elseif fun <= 25 
    err = 1e-1; 
elseif fun == 48 
    err = 1e-13; 
elseif fun == 49 
    err = 0.00001; 
elseif fun == 50 
    err = 0.0001; 
elseif fun == 51 
    err = 1e-10; 
else 
    err = 1e-4; 
end 
 
end
