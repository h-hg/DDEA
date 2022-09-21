% This MATLAB R2016b code is across minimization problems. 
% Please cite this article as: P.Liao, C.Sun, G.Zhang et al., Multi-surrogate multi-tasking optimization of expensive problems, 
% Knowledge-Based Systems (2020) 106262, https://doi.org/10.1016/j.knosys.2020.106262.
% For suggestions please contact: Peng Liao (Email: pengliao.lp@gmail.com)

%% The parameters of multi-tasking optimization are as follows:
pop_M=100; % population size 
gen=50; % generation count 
rmp=0.3; % random mating probability
%% The parameters of MS-MTO are as follows:
reps =30; % repetitions 30
dim=10;  % dimensionality

global initial_flag;
rand('twister',mod(floor(now*8640000),2^31-1));
addpath(genpath(pwd))
warning off 
for index =1:8       % F1:Ellipsoid  F2:Rosenbrock  F3:Ackley  F4:Griewank  F5:Rastrigin  F6:CEC05 F10  F7:CEC05 F16  F8:CEC05 F19                                                       
    initial_flag=0;
    Tasks = benchmark(index,dim);   
    fit=11*Tasks(1).dims;   % when D=10、20、30
%     fit=1000;             % when D=50、100、200
    achieve_num=Tasks(1).dims*2; % Initial size of archive(2*D)   
    data_MS_MTO(index)=MS_MTO(Tasks,pop_M,gen,rmp,reps,fit,achieve_num);

end

    save('MS_MTO.mat','data_MS_MTO');