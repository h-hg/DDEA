clc;
clear;
close;
addpath(genpath(pwd));
warning off all
% 'S_ELLIPSOID'; 'S_ROSENBROCK';  'S_ACKLEY'; 'S_GRIEWANK'; 'S_RASTRIGIN';  'CEC05_F10'; 'CEC05_F16';  'CEC05_F19';
funcArr={'S_ELLIPSOID'; 'S_ROSENBROCK';  'S_ACKLEY';  'S_GRIEWANK'; 'S_RASTRIGIN';  'CEC05_F10'; 'CEC05_F16';  'CEC05_F19';}; % functions
dimsArr = [10 30 50 100]; % dimension
runs = 20; % run times
o = length(funcArr); 
d = size(dimsArr,2);

for i = 1:d
    dims = dimsArr(i);
    for j = 1:o
        fname = cell2mat(funcArr(j));  
        result = [];
        time = [];
        for r = 1:runs
            % experiment parameters
            FUN=@(x) feval(fname,x); 
            [Xmin, Xmax] = variable_domain(fname); 
            LB = repmat((Xmin),1,dims);
            UB = repmat((Xmax),1,dims,1);
            initial_DS = 11*dims;
            hx=repmat(LB,initial_DS,1)+(repmat(UB,initial_DS,1)-repmat(LB,initial_DS,1)).*lhsdesign(initial_DS,dims);
            hf = FUN(hx)';
            DATA = [hx, hf]; % initial Data
            model_num = 4;  
            t1 = clock; 
            
            % prepare good fitness data and center point
            [DB_f,id]=sort(hf); 
            DB_x = hx(id,:);
            DB = [DB_x DB_f];
            DS = length(hf);
            GDS = floor(DS*0.2); % top 20% in the DB to form GD
            if GDS < 10 
                GDS = 10;
            end
            if GDS > DS 
                GDS = DS;
            end
            GD_x = DB_x(1:GDS,:);
            center = mean(GD_x,1);

            % train and test dataset
            E = zeros(1,model_num); % Model Error Criterion
            R = zeros(1,model_num); % Distance Deviation Criterion
            test_subset = 10; % number of subset
            for j = 1:test_subset
                % dataset divide
                index = randperm(GDS,floor(GDS/2));
                test_hf = DB_f(index); test_hx = DB_x(index,:);
                train_hf = DB_f;  train_hx = DB_x;
                train_hf(index) = [];  train_hx(index,:) = [];
                DS_train = length(train_hf);
                
                % RBF parameters
                ghxd=real(sqrt(train_hx.^2*ones(size(train_hx'))+ones(size(train_hx))*(train_hx').^2-2*train_hx*(train_hx')));
                spr=max(max(ghxd))/(dims*DS_train)^(1/dims);
                
                % train models
                for i = 1:model_num
                    % RBF network
                    h = 4*(i-1);
                    spr= spr * 2^h;
                    net=newrbe(train_hx',train_hf',spr);
                    RBF_FUN=@(x) sim(net,x');
                    % test error
                    error = abs(test_hf - RBF_FUN(test_hx)');
                    E(i) = E(i) + sum(error);
                end
            end
            E = E/test_subset;

            % build models including all samples for predict
            predict = zeros(1,model_num);
            predict_pos = [];
            for i = 1:model_num
                % RBF network use all data
                ghxd=real(sqrt(hx.^2*ones(size(hx'))+ones(size(hx))*(hx').^2-2*hx*(hx')));
                spr=max(max(ghxd))/(dims*DS)^(1/dims);
                h = 4*(i-1);
                spr= spr * 2^h;
                net=newrbe(hx',hf',spr);
                RBF_FUN=['RBF_FUN',num2str(i)];
                eval([RBF_FUN,'=@(x) sim(net,x'');']);
                eval(['RBF_FUN=',RBF_FUN,';']);

                % model predict
                maxgen = 30000+300*dims;
                minerror = 1e-10;
                RBF_FUN_Arr{i} = RBF_FUN;
                [best_pos,bestever] = SLPSO(dims, maxgen, RBF_FUN, minerror, GD_x);
                
                % obtain real fitness to verify effect of this algorithm
                predict_pos(i,:) = best_pos;
                predict(i) = FUN(best_pos);
            end
            
            
            T = 10;
            for k = 1:T
                R_temp = pdist2(DB_x(k,:), predict_pos);
                R = R_temp + R;
            end
            
            R_N = mapminmax(R);
            E_N = mapminmax(E);
            loss = (R_N + E_N)/2;
            [~, index]=min(loss);
            
            model_selected = index; % selected model index
            pos = predict_pos(index,:); % predicted solution
            pred = predict(index); % fitness of solution

            predict_N = mapminmax(predict);
            result(r) = pred    % result
            
            totaltime = etime(clock,t1);
            time(r) = totaltime;
        end
        best_result = min(result);
        worst_result = max(result);
        mean_result = mean(result);
        median_result = median(result);
        std_result    = std(result);
        mean_time = mean(time);
        out1 = [best_result,worst_result,mean_result,median_result,std_result];
        save(strcat('result/',fname,' runs=',num2str(runs),' dims=',num2str(dims)));
    end
end
