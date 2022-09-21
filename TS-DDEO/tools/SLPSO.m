%------------------------------------------------------------------------
% This code is part of the program that produces the results in the following paper:
% Huixiang Zhen, Wenyin Gong, Ling Wang, Fei Ming, and Zuowen Liao. "Two-stage Data-driven Evolutionary Optimization for High-dimensional Expensive Problems", IEEE Transactions on Cybernetics, accepted, 2021.
% You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%----------------------------------------------------------------------------------------------------------------------------------------

function [best_pos,bestever] = SLPSO(d, maxgen,gmd,minerror,ghx)
time_begin=tic;
% disp('SLPSO global search');
n = d;
M = 100; beta=0.01;
m = M + fix(d/10);
c3 = d/M*beta;
PL = zeros(m,1);

for i = 1 : m
    PL(i) = (1 - (i - 1)/m)^log(sqrt(ceil(d/M)));
end

%initialization
p = zeros(m, d); 
lu = [min(ghx); max(ghx)];

XRRmin = repmat(lu(1, :), m, 1);
XRRmax = repmat(lu(2, :), m, 1);
p = XRRmin + (XRRmax - XRRmin) .* lhsdesign(m, d);
v = zeros(m,d);
bestever = 1e200;
best_pos  = zeros(1,d);

FES = 0;
gen = 0;
flag_er=0;

tic;
%main loop
% while(FES < maxfe)
    while(gen < maxgen)  
        best_old=bestever;      
        bestpos_old=best_pos;   
        FES = FES + m;
        fitness = gmd(p);   
        %population sorting
        [fitness rank] = sort(fitness, 'descend'); 
        p = p(rank,:);
        v = v(rank,:);
        besty = fitness(m);
        bestp = p(m, :);
        [bestever,id] = min([besty, bestever]);
        best_new = bestever;    
        if id == 1
            best_pos = bestp;      
            bestpos_new=best_pos;  
        elseif id == 2
            best_pos = bestpos_old;
            bestpos_new=best_pos;
        end

        error=abs(best_old-best_new);
        if error <= minerror
            flag_er=flag_er+1;
        else
            flag_er=0;
        end
        if flag_er >=20
            break;
        end

        center = ones(m,1)*mean(p);
        randco1 = rand(m, d);
        randco2 = rand(m, d);
        randco3 = rand(m, d);
        winidxmask = repmat([1:m]', [1 d]);
        winidx = winidxmask + ceil(rand(m, d).*(m - winidxmask));
        pwin = p;
        for j = 1:d
            pwin(:,j) = p(winidx(:,j),j);
        end
        %learning
         lpmask = repmat(rand(m,1) < PL, [1 d]);
         lpmask(m,:) = 0;
         v1 =  1*(randco1.*v + randco2.*(pwin - p) + c3*randco3.*(center - p));
         p1 =  p + v1; 
         v = lpmask.*v1 + (~lpmask).*v;         
         p = lpmask.*p1 + (~lpmask).*p;
         %boundary
        for i = 1:m
            p(i,:) = max(p(i,:), lu(1,:));
            p(i,:) = min(p(i,:), lu(2,:));
        end
        gen = gen + 1;
    end
    time_cost=toc(time_begin);
%     disp(['local search time_cost=' num2str(time_cost)]);
end
