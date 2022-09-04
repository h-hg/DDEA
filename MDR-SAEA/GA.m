% Input:
% Parent         -Parent individuals for the evolutionary process
% lower          -Lower bound of decision variables 
% upper          -Upper bound of decision variables
%
% Output:
% Offspring           -Offspring individuals
%%%% Authors: Zheng Tan, Handing Wang, Yaochu Jin
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


function Offspring = GA(Parent,lower,upper)
    proC = 1;
    disC = 20;
    proM = 1;
    disM = 20;

    Parent1 = Parent(1:floor(end/2),:);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);

    %% Genetic operators for real encoding
    % Simulated binary crossover
    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                 (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
    % Polynomial mutation
    Lower = repmat(lower,2*N,1);
    Upper = repmat(upper,2*N,1);
    Site  = rand(2*N,D) < proM/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
    end

