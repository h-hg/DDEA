function Tasks = benchmark(index,dim)
%BENCHMARK function
%   Input
%   - index: the index number of problem set
%   - dim: the dimensions of the problem set
%   Output:
%   - Tasks: benchmark problem set
%   - task(1): global model
%   - task(2): local model
%   - task(3): Real function

    switch(index)
        
         case 1 
             
            Tasks(1).dims = dim;    % dimensionality of model 1
%           Tasks(1).fnc = @(x)model_1(x);
            Tasks(1).Lb=-5.12*ones(1,dim);   % Upper bound of model 1
            Tasks(1).Ub=5.12*ones(1,dim);    % Lower bound of model 1
            
            Tasks(2).dims = dim;    % dimensionality of model 2
%           Tasks(2).fnc = @(x)model_2(x);
            Tasks(2).Lb=-5.12*ones(1,dim);   % Upper bound of model 2
            Tasks(2).Ub=5.12*ones(1,dim);    % Lower bound of model 2
            
            Tasks(3).dims = dim;    % dimensionality of real function
            Tasks(3).fnc = @(x)Ellipsoid(x);
            Tasks(3).Lb=-5.12*ones(1,dim);   % Upper bound of real function
            Tasks(3).Ub=5.12*ones(1,dim);    % Lower bound of real function
                           
        case 2   

            Tasks(1).dims = dim;    % dimensionality of model 1
%           Tasks(1).fnc = @(x)model_1(x);
            Tasks(1).Lb=-2.048*ones(1,dim);   % Upper bound of model 1
            Tasks(1).Ub=2.048*ones(1,dim);    % Lower bound of model 1
            
            Tasks(2).dims = dim;    % dimensionality of model 2
%           Tasks(2).fnc = @(x)model_2(x);
            Tasks(2).Lb=-2.048*ones(1,dim);   % Upper bound of model 2       
            Tasks(2).Ub=2.048*ones(1,dim);    % Lower bound of model 2
            
            Tasks(3).dims = dim;    % dimensionality of real function
            Tasks(3).fnc = @(x)Rosenbrock(x);
            Tasks(3).Lb=-2.048*ones(1,dim);   % Upper bound of real function
            Tasks(3).Ub=2.048*ones(1,dim);    % Lower bound of real function
            
        case 3      

            Tasks(1).dims = dim;    % dimensionality of model 1
%           Tasks(1).fnc = @(x)model_1(x);
            Tasks(1).Lb=-32.768*ones(1,dim);   % Upper bound of model 1
            Tasks(1).Ub=32.768*ones(1,dim);    % Lower bound of model 1
            
            Tasks(2).dims = dim;    % dimensionality of model 2
%           Tasks(2).fnc = @(x)model_2(x);
            Tasks(2).Lb=-32.768*ones(1,dim);   % Upper bound of model 2
            Tasks(2).Ub=32.768*ones(1,dim);    % Lower bound of model 2
            
            Tasks(3).dims = dim;    % dimensionality of real function
            Tasks(3).fnc = @(x)Ackley(x);
            Tasks(3).Lb=-32.768*ones(1,dim);   % Upper bound of real function
            Tasks(3).Ub=32.768*ones(1,dim);    % Lower bound of real function
               
       case 4
         
            Tasks(1).dims = dim;    % dimensionality of model 1
%           Tasks(1).fnc = @(x)model_1(x);
            Tasks(1).Lb=-600*ones(1,dim);   % Upper bound of model 1
            Tasks(1).Ub=600*ones(1,dim);    % Lower bound of model 1
            
            Tasks(2).dims = dim;    % dimensionality of model 2
%           Tasks(2).fnc = @(x)model_2(x);
            Tasks(2).Lb=-600*ones(1,dim);   % Upper bound of model 2
            Tasks(2).Ub=600*ones(1,dim);    % Lower bound of model 2
            
            Tasks(3).dims = dim;    % dimensionality of real function
            Tasks(3).fnc = @(x)Griewank(x);
            Tasks(3).Lb=-600*ones(1,dim);   % Upper bound of real function
            Tasks(3).Ub=600*ones(1,dim);    % Lower bound of real function
                   
        case 5
           
            Tasks(1).dims = dim;    % dimensionality of model 1
%           Tasks(1).fnc = @(x)model_1(x);
            Tasks(1).Lb=-5.12*ones(1,dim);   % Upper bound of model 1
            Tasks(1).Ub=5.12*ones(1,dim);    % Lower bound of model 1
            
            Tasks(2).dims = dim;    % dimensionality of model 2
%           Tasks(2).fnc = @(x)model_2(x);
            Tasks(2).Lb=-5.12*ones(1,dim);   % Upper bound of model 2
            Tasks(2).Ub=5.12*ones(1,dim);    % Lower bound of model 2
             
            Tasks(3).dims = dim;    % dimensionality of real function
            Tasks(3).fnc = @(x)Rastrigin(x);
            Tasks(3).Lb=-5.12*ones(1,dim);   % Upper bound of real function
            Tasks(3).Ub=5.12*ones(1,dim);    % Lower bound of real function
           	
        case 6
          
            Tasks(1).dims = dim;    % dimensionality of model 1
%           Tasks(1).fnc = @(x)model_1(x);
            Tasks(1).Lb=-5*ones(1,dim);   % Upper bound of model 1
            Tasks(1).Ub=5*ones(1,dim);    % Lower bound of model 1
            
            Tasks(2).dims = dim;    % dimensionality of model 2
%           Tasks(2).fnc = @(x)model_2(x);
            Tasks(2).Lb=-5*ones(1,dim);   % Upper bound of model 2
            Tasks(2).Ub=5*ones(1,dim);    % Lower bound of model 2
            
            Tasks(3).dims = dim;    % dimensionality of real function
            Tasks(3).fnc = @(x)benchmark_func(x,10);
            Tasks(3).Lb=-5*ones(1,dim);   % Upper bound of real function
            Tasks(3).Ub=5*ones(1,dim);    % Lower bound of real function
     
       case 7 

            Tasks(1).dims = dim;    % dimensionality of model 1
%           Tasks(1).fnc = @(x)model_1(x);
            Tasks(1).Lb=-5*ones(1,dim);   % Upper bound of model 1
            Tasks(1).Ub=5*ones(1,dim);    % Lower bound of model 1
            
            Tasks(2).dims = dim;    % dimensionality of model 2
%           Tasks(2).fnc = @(x)model_2(x);
            Tasks(2).Lb=-5*ones(1,dim);   % Upper bound of model 2
            Tasks(2).Ub=5*ones(1,dim);    % Lower bound of model 2
             
            Tasks(3).dims = dim;    % dimensionality of real function
            Tasks(3).fnc = @(x)benchmark_func(x,16);
            Tasks(3).Lb=-5*ones(1,dim);   % Upper bound of real function
            Tasks(3).Ub=5*ones(1,dim);    % Lower bound of real function
               
       case 8 

            Tasks(1).dims = dim;    % dimensionality of model 1
%           Tasks(1).fnc = @(x)model_1(x);
            Tasks(1).Lb=-5*ones(1,dim);   % Upper bound of model 1
            Tasks(1).Ub=5*ones(1,dim);    % Lower bound of model 1
            
            Tasks(2).dims = dim;    % dimensionality of model 2
%           Tasks(2).fnc = @(x)model_2(x);
            Tasks(2).Lb=-5*ones(1,dim);   % Upper bound of model 2
            Tasks(2).Ub=5*ones(1,dim);    % Lower bound of model 2
             
            Tasks(3).dims = dim;    % dimensionality of real function
            Tasks(3).fnc = @(x)benchmark_func(x,19);
            Tasks(3).Lb=-5*ones(1,dim);   % Upper bound of real function
            Tasks(3).Ub=5*ones(1,dim);    % Lower bound of real function
                    
    end
end