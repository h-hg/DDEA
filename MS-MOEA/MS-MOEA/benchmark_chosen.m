% An Adaptive Model Switch-based Surrogate-Assisted Evolutionary Algorithm 
% for Noisy Expensive Multi-Objective Optimization
%------------------------------- Reference --------------------------------
% N. Zheng, H. Wang, and B. Yuan, An Adaptive Model Switch-based 
% Surrogate-Assisted Evolutionary Algorithm for Noisy Expensive 
% Multi-Objective Optimization, Complex & Intelligent Systems, accepted, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 HandingWangXD Group. Permission is granted to copy and
% use this code for research, noncommercial purposes, provided this
% copyright notice is retained and the origin of the code is cited. The
% code is provided "as is" and without any warranties, express or implied.

% This code is written by Nan Zheng.
% Email: nanszheng@163.com
%选择测试问题
function [num_vari,num_obj,Size,max_evaluation,decision_space,PF]=benchmark_chosen(fun_test)
switch fun_test
     case 'DTLZ1M2'
    
        num_vari=7;%决策变量维数
        num_obj=2;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造     
        Q=load('DTLZ1M2.mat');
        PF=Q.P;
     case 'DTLZ1M3'
    
        num_vari=7;%决策变量维数
        num_obj=3;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造     
        Q=load('DTLZ1M3.mat');
        PF=table2array(Q.DTLZ1);        
     case 'DTLZ2M2'
    
        num_vari=12;%决策变量维数
        num_obj=2;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造      
        Q=load('DTLZ2M2.mat');
        PF=Q.P;
     case 'DTLZ2M3'
    
        num_vari=12;%决策变量维数
        num_obj=3;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造      
        Q=load('DTLZ2M3.mat');
        PF=table2array(Q.DTLZ2);        
     case 'DTLZ3M2'
    
        num_vari=7;%决策变量维数
        num_obj=2;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造
        Q=load('DTLZ3M2.mat');
        PF=Q.P;
     case 'DTLZ3M3'
    
        num_vari=7;%决策变量维数
        num_obj=3;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造
        Q=load('DTLZ3M3.mat');
        PF=table2array(Q.DTLZ3);        
     case 'DTLZ4M2'
    
        num_vari=12;%决策变量维数
        num_obj=2;%目标数
        Size=100;
        max_evaluation=300;       
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造 
        Q=load('DTLZ4M2.mat');
        PF=Q.P;
     case 'DTLZ4M3'
    
        num_vari=12;%决策变量维数
        num_obj=3;%目标数
        Size=100;
        max_evaluation=300;       
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造 
        Q=load('DTLZ4M3.mat');
        PF=table2array(Q.DTLZ4);        
      case 'DTLZ5M2'
    
        num_vari=12;%决策变量维数
        num_obj=2;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造   
        Q=load('DTLZ5M2.mat');
        PF=Q.P;
      case 'DTLZ5M3'
    
        num_vari=12;%决策变量维数
        num_obj=3;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造   
        Q=load('DTLZ5M3.mat');
        PF=table2array(Q.DTLZ5);        
      case 'DTLZ6M2'
    
        num_vari=12;%决策变量维数
        num_obj=2;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造   
        Q=load('DTLZ6M2.mat');
        PF=Q.P;
      case 'DTLZ6M3'
    
        num_vari=12;%决策变量维数
        num_obj=3;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造   
        Q=load('DTLZ6M3.mat');
        PF=table2array(Q.DTLZ6);        
      case 'DTLZ7M2'
    
        num_vari=12;%决策变量维数
        num_obj=2;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造      
        Q=load('DTLZ7M2.mat');
        PF=Q.P;
      case 'DTLZ7M3'
    
        num_vari=12;%决策变量维数
        num_obj=3;%目标数
        Size=100;
        max_evaluation=300;        
        decision_space=[zeros(1,num_vari);ones(1,num_vari)];%决策空间构造      
        Q=load('DTLZ7M3.mat');
        PF=table2array(Q.DTLZ7);        
          
        
    otherwise
        error('objective function is not defined!');        
end   
end