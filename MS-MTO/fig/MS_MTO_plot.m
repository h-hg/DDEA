% This MATLAB R2016b code is across minimization problems. 
% This is the data of MS-MTO running in this paper.
% Test set is CEC2005 and CEC2017.
% F1:Ellipsoid  F2:Rosenbrock  F3:Ackley  F4:Griewank  F5:Rastrigin  F6:CEC2005 F10  F7:CEC2005 F16  F8:CEC2005 F19 
% F9:CEC2017 F1  F10:CEC2017 F4  F11:CEC2017 F11  F12:CEC2017 F21
% Each function runs independently 30 times
% When dim=20 or 200, index=1,2,3,4,5,6,7,8
% When dim=10, 30, 50 or 100, index=1,2,3,4,5,6,7,8,9,10,11,12
% Data format description:
% -wall_clock_time:Single run time
% -best_one:The final output solution after running the algorithm
% -Toall_FES:The number of times the real function is called in the algorithm run
% -Toall_BestFitss:The best fitness of each generation in the running of the algorithm
% For suggestions please contact: Peng Liao (Email: pengliao.lp@gmail.com)
clear all

index=1;
dim=100; 
runs=30;

if dim==10
    load('MS_MTO_10.mat')
    x_=110;
    x_10=[20:10:x_];%MS-MTO
    for   run=1:runs          
        count=1;  
        for j=1:size(MS_MTO_10(index).Toall_BestFitss(run,:),2)
           while count<=MS_MTO_10(index).FES(1,j)  
           data_10(run,count)=MS_MTO_10(index).Toall_BestFitss(run,j);
           count=count+1;
           end
        end
    end
    data_10(34,:)=mean(data_10(1:runs,:));
    data_10(35,:)=std(data_10(1:runs,:));
    figure(1)
    plot(x_10,((data_10(34,x_10(1,:)))),'-r','linewidth',2);
    
elseif dim==20
    load('MS_MTO_20.mat')
    x_=220;
    x_20=[40:20:x_];%MS-MTO
    for   run=1:runs          
        count=1;  
        for j=1:size(MS_MTO_20(index).Toall_BestFitss(run,:),2)
           while count<=MS_MTO_20(index).FES(1,j)  
           data_20(run,count)=MS_MTO_20(index).Toall_BestFitss(run,j);
           count=count+1;
           end
        end
    end
    data_20(34,:)=mean(data_20(1:runs,:));
    data_20(35,:)=std(data_20(1:runs,:));
    figure(1)
    plot(x_20,((data_20(34,x_20(1,:)))),'-r','linewidth',2);
    
elseif dim==30
    load('MS_MTO_30.mat')
    x_=330;
    x_30=[60:30:x_];%MS-MTO
    for   run=1:runs          
        count=1;  
        for j=1:size(MS_MTO_30(index).Toall_BestFitss(run,:),2)
           while count<=MS_MTO_30(index).FES(1,j)  
           data_30(run,count)=MS_MTO_30(index).Toall_BestFitss(run,j);
           count=count+1;
           end
        end
    end
    data_30(34,:)=mean(data_30(1:runs,:));
    data_30(35,:)=std(data_30(1:runs,:));
    figure(1)
    plot(x_30,((data_30(34,x_30(1,:)))),'-r','linewidth',2);
    
elseif dim==50
    load('MS_MTO_50.mat')
    x_=1000;
    x_50=[100:50:x_];%MS-MTO
    for   run=1:runs          
        count=1;  
        for j=1:size(MS_MTO_50(index).Toall_BestFitss(run,:),2)
           while count<=MS_MTO_50(index).FES(1,j)  
           data_50(run,count)=MS_MTO_50(index).Toall_BestFitss(run,j);
           count=count+1;
           end
        end
    end
    data_50(34,:)=mean(data_50(1:runs,:));
    data_50(35,:)=std(data_50(1:runs,:));
    figure(1)
    plot(x_50,((data_50(34,x_50(1,:)))),'-r','linewidth',2);
    
elseif dim==100
    load('MS_MTO_100.mat')
    x_=1000;
    x_100=[200:50:x_];%MS-MTO
    for   run=1:runs          
        count=1;  
        for j=1:size(MS_MTO_100(index).Toall_BestFitss(run,:),2)
           while count<=MS_MTO_100(index).FES(1,j)  
           data_100(run,count)=MS_MTO_100(index).Toall_BestFitss(run,j);
           count=count+1;
           end
        end
    end
    data_100(34,:)=mean(data_100(1:runs,:));
    data_100(35,:)=std(data_100(1:runs,:));
    figure(1)
    plot(x_100,((data_100(34,x_100(1,:)))),'-r','linewidth',2);
    
elseif dim==200
    load('MS_MTO_200.mat')
    x_=1000;
    x_200=[400:50:x_];%MS-MTO
    for   run=1:runs          
        count=1;  
        for j=1:size(MS_MTO_200(index).Toall_BestFitss(run,:),2)
           while count<=MS_MTO_200(index).FES(1,j)  
           data_200(run,count)=MS_MTO_200(index).Toall_BestFitss(run,j);
           count=count+1;
           end
        end
    end
    data_200(34,:)=mean(data_200(1:runs,:));
    data_200(35,:)=std(data_200(1:runs,:));
    figure(1)
    plot(x_200,((data_200(34,x_200(1,:)))),'-r','linewidth',2);
    
end


  
    xlim([0 x_]);
    xlabel('Number of fitness evaluation') ;
%     ylabel('Mean fitness value');
    ylabel('Mean fitness value (Natural Log)');
  
    h= legend('MS-MTO');
    set(h,'Fontsize',18);
    set(gca,'box','on','Fontsize',18);

  
    
    
    



