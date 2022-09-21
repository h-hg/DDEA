clc
clear all
cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));
G_Optimal = [];
Iter_OptimalResult = [];
Time = [];
run_times = 20;

max_evaluations = 1000;

for i=5:21
    if i<=8 || i==21
        [xmin xmax vmin vmax] = bound_definition(i);
        if (i>=1 && i<= 4) || (i>=9 && i<=10)
            dimension = 50;
        else if (i>=5 && i<= 8) || (i>=11 && i<=12) || i==21
                dimension = 30;
            else
                dimension = 1000;
            end
        end
        popsize = 50;
        %[G_Optimal1 Iter_OptimalResult1 Time1] = SPSO_original(i,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations);
        %[G_Optimal1 Iter_OptimalResult1 Time1] = SPSO_version2(i,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations);
        %[G_Optimal1 Iter_OptimalResult1 Time1] = SPSO_w(i,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations);
        [G_Optimal1 Iter_OptimalResult1 Time1] = SPSO_w(i,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations);
        %[G_Optimal1 Iter_OptimalResult1 Time1] = oSPSO_w(i,xmin,xmax,vmin,vmax,run_times,dimension,popsize,max_evaluations);
        for r=1:run_times
            G_Optimal(i,r,:) = G_Optimal1(r,:);
                   Iter_OptimalResult(i,r,:) = Iter_OptimalResult1(r,:);
            Time(i,r,:) = Time1(r,:);
        end
    end
    if i==5
        for r=1:run_times
            iter_result1(r,:) = Iter_OptimalResult1(r,:);
        end
    else if i == 6
            for r=1:run_times
                iter_result2(r,:) = Iter_OptimalResult1(r,:);
            end
        else if i==7
                for r=1:run_times
                    iter_result3(r,:) = Iter_OptimalResult1(r,:);
                end
            else if i==8
                    for r=1:run_times
                        iter_result4(r,:) = Iter_OptimalResult1(r,:);
                    end
                else if i==21
                        for r=1:run_times
                            iter_result5(r,:) = Iter_OptimalResult1(r,:);
                        end
                    end
                end
            end
        end
    end
end

