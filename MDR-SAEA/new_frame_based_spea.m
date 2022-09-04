%%%% Authors: Zheng Tan, Handing Wang, Shulei Liu
%%%% Xidian University, China and Chinese Academy of Military Science, China.
%%%% EMAIL: zhengtan@stu.xidian.edu.cn, hdwang @ xidian.edu.cn
%%%% WEBSITE: https://sites.google.com/site/handingwanghomepage
%%%% DATE: March 2021
% ------------------------------------------------------------------------
% This code is part of the program that produces the results in the following paper:
%
% Zheng Tan, Handing Wang, Shulei Liu, Multi-Stage Dimension Reduction for Expensive Sparse Multi-Objective Optimization Problems, Neurocomputing, vol.440, no.14, pp.159–174, 2021.
%
% You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
% ------------------------------------------------------------------------

%% Input
    global evaluation; % Number of expensive function evaluations used
    evaluation = 0 ;
    M = 2; % Number of objective
    D = 1000; % Dimension of decision variables
    
    % Lower and upper bound for decision variables
    lower    = [zeros(1,M-1)+0,zeros(1,D-M+1)-1];
    upper    = [zeros(1,M-1)+1,zeros(1,D-M+1)+2];
    
    % Objective function handle 
    g_1 = @g1; 
    f = @CalObj;
    g_2 = @g2;
    
    sparsity = 0.01; %sparsity of the problem
    
    %% NDS-based Dimension Selection
    [dim,TMask,TempPop,popobj] = dim_selection(f,sparsity,g_1,g_2,D,lower,upper);
%% Initialization
   
    N = size(dim,2);
    fitness  = zeros(1,N);
    for i = 1:3
            fitness    = fitness + NDSort(popobj((i-1)*N+1:i*N,:),inf);
    end
    init_time = 3;
    if init_time*size(dim,2)>0.5*D
        init_time = floor(0.5*D/size(dim,2));
    end
    for i = 1 :init_time
        Dec = lower+(upper-lower).*rand(N,D);
        Mask_temp       = eye(N);
        Mask = zeros(N,D);
        for j =1:N
            Mask(j,dim(find(Mask_temp(j,:)))) = 1;
        end
        population = Dec.*Mask;
        PopObj     = f(population,sparsity,g_1,g_2);
        popobj     = [popobj;PopObj];
        TMask      = [TMask;Mask_temp];
        TempPop    = [TempPop;population];
        fitness    = fitness + NDSort(PopObj,inf);
    end
    % Generate initial population
    Dec = lower+(upper-lower).*rand(N,D);
    Mask_dim = eye(N);
    Mask = zeros(N,D);
    for i =1:N
        Mask(i,dim(find(Mask_dim(i,:)))) = 1;
    end
    Population = Dec.*Mask;
    PopObj_1     = f(Population,sparsity,g_1,g_2);
    PopObj = [popobj;PopObj_1];
    num_e = 100;
    [population,mask_dim,FrontNo,CrowdDis,popobj] = EnvironmentalSelection([TempPop;Population],PopObj,[TMask;Mask_dim],num_e);
    
    
%% Search Proper Dimension

num_time = 0;
non_zero_all = zeros(1,50);
non_zero = zeros(50,num_e);
for time = 1:50
    MatingPool       = TournamentSelection(2,2*100,FrontNo,-CrowdDis);
    mask_off = operate_spea(mask_dim(MatingPool,:),1,fitness);
    Dec = lower+(upper-lower).*rand(num_e,D);
    Mask = zeros(num_e,D);
    for i =1:num_e
        Mask(i,dim(find(mask_off(i,:)))) = 1;
    end
    population_off = Dec.*Mask;
    PopObj = [popobj;f(population_off,sparsity,g_1,g_2)];
    [population,mask_dim,FrontNo,CrowdDis,popobj] = EnvironmentalSelection([population;population_off],PopObj,[mask_dim;mask_off],100);
    for i =1:num_e
        non_zero_temp = find(mask_dim(i,:));
        non_zero(time,i) = size(non_zero_temp,2);
    end
    num_temp = mode(non_zero(time,:),2);
    if time>1
    if num_temp == dim_init
        num_time = num_time+1;
    else
        num_time = 1;
    end
    end
    dim_init = num_temp;
    if num_time>4
        break
    end
end

%% Selection of Dimension For SAEA

dim_saea = dim_base(mask_dim,dim,num_temp); 
result = saea_check(dim_saea,num_temp);
if result>0
    dim_saea = dim_final_selection(dim_saea,f,sparsity,g_1,g_2,num_temp);
end

%% Output
% num_temp: Obtained number of non-zero dimensions
% dim_saea: Obtained non-zero dimensions



        %% Calculate objective values

% SMOP1 benchmark problem
function PopObj = CalObj(X,sparsity,g_1,g_2)
    global evaluation;
    [N,D] = size(X);
    evaluation = evaluation + N;
    M = 2;
    K = ceil(sparsity*(D-M+1));
    g = sum(g_1(X(:,M:M+K-1),pi/3),2) + sum(g_2(X(:,M+K:end),0),2);
    PopObj = repmat(1+g/(D-M+1),1,M).*fliplr(cumprod([ones(N,1),X(:,1:M-1)],2)).*[ones(N,1),1-X(:,M-1:-1:1)];
end



function g = g1(x,t)
    g = (x-t).^2;
end

function g = g2(x,t)
    g = 2*(x-t).^2 + sin(2*pi*(x-t)).^2;
end

function g = g3(x,t)
    g = 4-(x-t)-4./exp(100*(x-t).^2);
end

function g = g4(x,t)
    g = (x-pi/3).^2 + t.*sin(6*pi*(x-pi/3)).^2;
end


