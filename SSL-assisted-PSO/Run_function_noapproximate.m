clc
clear all
G_Optimal = [];
Iter_OptimalResult = [];
Time = [];
run_times = 20;
popsize = 60;
max_evaluations = 1000;

for i=10:10
    [xmin xmax vmin vmax] = bound_definition(i);
    if (i>=1 && i<= 4) || (i>=9 && i<=10)
        dimension = 50;
    else if (i>=5 && i<= 8) || (i>=11 && i<=12)
            dimension = 100;
        else
            dimension = 1000;
        end
    end
    [G_Optimal1 Iter_OptimalResult1 Time1] = SPSO_noapproximate(i,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations);
    for r=1:run_times
        G_Optimal(i,r,:) = G_Optimal1(r,:);
        Iter_OptimalResult(i,r,:) = Iter_OptimalResult1(r,:);
        Time(i,r,:) = Time1(r,:);
    end
end
% 
% %%%global best
% for r=1:20
%     G_Optimal11(r,:) = G_Optimal(1,r,:);
%     G_Optimal22(r,:) = G_Optimal(2,r,:);
%     G_Optimal33(r,:) = G_Optimal(3,r,:);
%     G_Optimal44(r,:) = G_Optimal(4,r,:);
%     G_Optimal55(r,:) = G_Optimal(5,r,:);
%     G_Optimal66(r,:) = G_Optimal(6,r,:);
%     G_Optimal77(r,:) = G_Optimal(7,r,:);
%     G_Optimal88(r,:) = G_Optimal(8,r,:);
%     G_Optimal99(r,:) = G_Optimal(9,r,:);
%     G_Optimal1010(r,:) = G_Optimal(10,r,:);
%     G_Optimal1111(r,:) = G_Optimal(11,r,:);
%     G_Optimal1212(r,:) = G_Optimal(12,r,:);
%     G_Optimal1313(r,:) = G_Optimal(13,r,:);
%     G_Optimal1414(r,:) = G_Optimal(14,r,:);
% end
% save 'e:\FI\results_noapproximate\G_Optimal11.txt' G_Optimal11 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal22.txt' G_Optimal22 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal33.txt' G_Optimal33 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal44.txt' G_Optimal44 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal55.txt' G_Optimal55 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal66.txt' G_Optimal66 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal77.txt' G_Optimal77 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal88.txt' G_Optimal88 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal99.txt' G_Optimal99 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal1010.txt' G_Optimal1010 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal1111.txt' G_Optimal1111 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal1212.txt' G_Optimal1212 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal1313.txt' G_Optimal1313 -ascii
% save 'e:\FI\results_noapproximate\G_Optimal1414.txt' G_Optimal1414 -ascii
% 
% %%%%%%% Best for each evaluation
% for r=1:20
%     Iter_OptimalResult11(r,:) = Iter_OptimalResult(1,r,:);
%     Iter_OptimalResult22(r,:) = Iter_OptimalResult(2,r,:);
%     Iter_OptimalResult33(r,:) = Iter_OptimalResult(3,r,:);
%     Iter_OptimalResult44(r,:) = Iter_OptimalResult(4,r,:);
%     Iter_OptimalResult55(r,:) = Iter_OptimalResult(5,r,:);
%     Iter_OptimalResult66(r,:) = Iter_OptimalResult(6,r,:);
%     Iter_OptimalResult77(r,:) = Iter_OptimalResult(7,r,:);
%     Iter_OptimalResult88(r,:) = Iter_OptimalResult(8,r,:);
%     Iter_OptimalResult99(r,:) = Iter_OptimalResult(9,r,:);
%     Iter_OptimalResult1010(r,:) = Iter_OptimalResult(10,r,:);
%     Iter_OptimalResult1111(r,:) = Iter_OptimalResult(11,r,:);
%     Iter_OptimalResult1212(r,:) = Iter_OptimalResult(12,r,:);
%     Iter_OptimalResult1313(r,:) = Iter_OptimalResult(13,r,:);
%     Iter_OptimalResult1414(r,:) = Iter_OptimalResult(14,r,:);
% end
% save 'e:\FI\results_noapproximate\F11.txt' Iter_OptimalResult11 -ascii
% save 'e:\FI\results_noapproximate\F22.txt' Iter_OptimalResult22 -ascii
% save 'e:\FI\results_noapproximate\F33.txt' Iter_OptimalResult33 -ascii
% save 'e:\FI\results_noapproximate\F44.txt' Iter_OptimalResult44 -ascii
% save 'e:\FI\results_noapproximate\F55.txt' Iter_OptimalResult55 -ascii
% save 'e:\FI\results_noapproximate\F66.txt' Iter_OptimalResult66 -ascii
% save 'e:\FI\results_noapproximate\F77.txt' Iter_OptimalResult77 -ascii
% save 'e:\FI\results_noapproximate\F88.txt' Iter_OptimalResult88 -ascii
% save 'e:\FI\results_noapproximate\F99.txt' Iter_OptimalResult99 -ascii
% save 'e:\FI\results_noapproximate\F1010.txt' Iter_OptimalResult1010 -ascii
% save 'e:\FI\results_noapproximate\F1111.txt' Iter_OptimalResult1111 -ascii
% save 'e:\FI\results_noapproximate\F1212.txt' Iter_OptimalResult1212 -ascii
% save 'e:\FI\results_noapproximate\F1313.txt' Iter_OptimalResult1313 -ascii
% save 'e:\FI\results_noapproximate\F1414.txt' Iter_OptimalResult1414 -ascii
% 
% %%%% Time
% for r=1:20
%     time11(r,:) = Time(1,r,:);
%     time22(r,:) = Time(2,r,:);
%     time33(r,:) = Time(3,r,:);
%     time44(r,:) = Time(4,r,:);
%     time55(r,:) = Time(5,r,:);
%     time66(r,:) = Time(6,r,:);
%     time77(r,:) = Time(7,r,:);
%     time88(r,:) = Time(8,r,:);
%     time99(r,:) = Time(9,r,:);
%     time1010(r,:) = Time(10,r,:);
%     time1111(r,:) = Time(11,r,:);
%     time1212(r,:) = Time(12,r,:);
%     time1313(r,:) = Time(13,r,:);
%     time1414(r,:) = Time(14,r,:);
% end
% save 'e:\FI\results_noapproximate\time11.txt' time11 -ascii
% save 'e:\FI\results_noapproximate\time22.txt' time22 -ascii
% save 'e:\FI\results_noapproximate\time33.txt' time33 -ascii
% save 'e:\FI\results_noapproximate\time44.txt' time44 -ascii
% save 'e:\FI\results_noapproximate\time55.txt' time55 -ascii
% save 'e:\FI\results_noapproximate\time66.txt' time66 -ascii
% save 'e:\FI\results_noapproximate\time77.txt' time77 -ascii
% save 'e:\FI\results_noapproximate\time88.txt' time88 -ascii
% save 'e:\FI\results_noapproximate\time99.txt' time99 -ascii
% save 'e:\FI\results_noapproximate\time1010.txt' time1010 -ascii
% save 'e:\FI\results_noapproximate\time1111.txt' time1111 -ascii
% save 'e:\FI\results_noapproximate\time1212.txt' time1212 -ascii
% save 'e:\FI\results_noapproximate\time1313.txt' time1313 -ascii
% save 'e:\FI\results_noapproximate\time1414.txt' time1414 -ascii
