% This is K-RVEA 
function [off,current_pop,Empty_ref_old] = evolve_K_RVEA(dmodel,Population,FunctionValue,Boundary,Int_Col,Vs,Empty_ref_old)

format compact;tic;
maxFE = 10000;
Coding = 'Real';
V = Vs;
N = size(V,1);
M = size(FunctionValue,2);
Fix_V = Vs;
up_var = M;
alpha = 1; 
cosineVV = V*V';
scosineVV = sort(cosineVV, 2, 'descend');
acosVV = acos(scosineVV(:,2));
refV = (acosVV);
% Empty_ref_old = 0;
FE = 0;
w=0;
ERR = zeros(size(Population,1),1);
metagen = round(maxFE/N) + 200;
while FE <=maxFE
%% SBX and polynomial mutation    
        MatingPool= F_mating(Population,N);
        Offspring = P_generator(MatingPool,Boundary,Coding,N);  
        Fitness = zeros(size(Offspring,1),M);
        ER = Fitness;
        for i = 1:size(Offspring,1)
            for k = 1:M
                [Fitness(i,k),~,ER(i,k)] = predictor(Offspring(i,:),dmodel(k));
            end
            
        end
        ER = sqrt(ER);
        FE = FE + size(Fitness,1);
        ER = mean(ER,2);
%         er_max_new= max(mean(er,2));
        Population = [Population; Offspring];       
        FunctionValue = [FunctionValue; Fitness];
        ERR = [ERR;ER];
        if M==2
            theta =  (w/(metagen))^alpha;
        else
             theta =  M*(w/(metagen))^alpha;
        end
        Selection= F_select(FunctionValue,V, theta, refV);
        Population = Population(Selection,:);
        FunctionValue = FunctionValue(Selection,:);
        ERR = ERR(Selection,:);
        
        if(mod(w, ceil(metagen*0.3)) == 0 && w > 0)
        
            Zmin = min(FunctionValue,[],1);	
            Zmax = max(FunctionValue,[],1);	

            V = Vs;

            V = V.*repmat((Zmax - Zmin)*1.0,N,1);
            tV = V;


            for i = 1:N
                V(i,:) = V(i,:)./norm(V(i,:));
            end;

            cosineVV = V*V';
            [scosineVV, neighbor] = sort(cosineVV, 2, 'descend');
            acosVV = acos(scosineVV(:,2));

            refV = (acosVV);

        end;
        w = w+1;
end; 
info_update = struct('c',{FunctionValue,V,theta,Fix_V,Empty_ref_old,refV,up_var,N,Population,ERR,Boundary,dmodel});
[off,Empty_ref_old] = update_metamodel(info_update);       
off = unique(off,'rows');
current_pop = Population;
end
                
              

