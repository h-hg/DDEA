% This is K-RVEA 
function [off,Empty_ref_old] = evolve_K_RVEA(model_ex,model_nex,Archive,Boundary,Empty_ref_old,id_ex,id_nex)
no_var = size(Boundary,2);
Population = Archive(:,1:no_var);
FunctionValue = Archive(:,no_var+1:end);

maxFE = 20000;
Coding = 'Real';

Vs = ref_vectors(2); % Generation of reference vectors
V = Vs;

N = size(V,1);
M = size(FunctionValue,2);
Fix_V = Vs;
up_var = M;
alpha = 0.5; 
cosineVV = V*V';
scosineVV = sort(cosineVV, 2, 'descend');
acosVV = acos(scosineVV(:,2));
refV = (acosVV);

FE = 0;
w=0;
ERR = zeros(size(Population,1),1);
% metagen = round(maxFE/N) + 200;
while FE <maxFE
%% SBX and polynomial mutation    
        MatingPool= F_mating(Population,N);
        Offspring = P_generator(MatingPool,Boundary,Coding,N); 
        
        ER = zeros(size(Offspring,1),M);
        
        Fitness = zeros(size(Offspring,1),M);
        
        [Fitness(:,id_ex),ER(:,id_ex)] = predict(model_ex,Offspring);
        
        [Fitness(:,id_nex),ER(:,id_nex)] = predict(model_nex,Offspring);
        
        FE = FE + size(Fitness,1);
        ER = mean(ER,2);

        Population = [Population; Offspring];       
        FunctionValue = [FunctionValue; Fitness];
        ERR = [ERR;ER];
        if M==2
            theta =  (FE/(maxFE))^alpha;
        else
             theta =  M*(w/(metagen))^alpha;
        end
        Selection= F_select(FunctionValue,V, theta, refV);
        Population = Population(Selection,:);
        FunctionValue = FunctionValue(Selection,:);
        ERR = ERR(Selection,:);
        
        if(mod(FE, ceil(maxFE*0.3)) == 0 && FE > 0)
        
            Zmin = min(FunctionValue,[],1);	
            Zmax = max(FunctionValue,[],1);	

            V = Vs;

            V = V.*repmat((Zmax - Zmin)*1.0,N,1);
            tV = V;


            for i = 1:N
                V(i,:) = V(i,:)./norm(V(i,:));
            end

            cosineVV = V*V';
            [scosineVV, neighbor] = sort(cosineVV, 2, 'descend');
            acosVV = acos(scosineVV(:,2));

            refV = (acosVV);

        end
        w = w+1;
end
info_update = struct('c',{FunctionValue,V,theta,Fix_V,Empty_ref_old,refV,up_var,N,Population,ERR,Boundary});
[off,Empty_ref_old] = update_metamodel(info_update);       
off = unique(off,'rows');
current_pop = Population;


P_archive = Archive(:,1:no_var);

Lia = ismember(off,P_archive,'rows');
r_unique = find(Lia(:,1)==0);
off = off(r_unique,:);
tt = 1;
while (isempty(off))
    rt = randperm(size(current_pop,1));
    off = current_pop (rt(1:1),:);
    Lia = ismember(off,P_archive,'rows');
    r_unique = find(Lia(:,1)==0);
    off = off(r_unique,:);
%                     size (off)
    tt = tt+1;
    if tt>1000
        break;
    end
end 
end
                
              

