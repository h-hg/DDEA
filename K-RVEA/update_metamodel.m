function [my_pop,Empty_ref_new] = update_metamodel(info)

FunctionValue = info(1).c;
V = info(2).c;
theta = info(3).c;
Fix_V = info(4).c;
Empty_ref_old = info(5).c;
refV = info(6).c;
up_var = info(7).c;
N = info(8).c;
Population = info(9).c;
ERR = info(10).c;

M = size(FunctionValue,2);
[~,Empty_ref_new] = reference (FunctionValue,Fix_V,refV,ERR,Population,theta);
delta = size(Empty_ref_new,1) - size(Empty_ref_old,1);

% Find the active reference vectors
[~,~,Non_empty_ref,APD_class,ER_class,Population_class] = reference (FunctionValue,V,refV,ERR,Population,theta);
% [~, ~, ~, ~ ,~,~,~,Empty_ref_new] = F_select(FunctionValue,V, theta, Fix_V,delta);

if size (Non_empty_ref,1) < up_var
    cluster_size = size(Non_empty_ref,1);
else 
    cluster_size = up_var;
end
% delta = size(Empty_ref_new,1) - size(Empty_ref_old,1);
if size(Population,1) <2
    my_pop = Population;
else
    [idx,temp_plot] = kmeans(Non_empty_ref, cluster_size,'start','uniform', 'emptyaction','singleton');
    my_pop = [];
    % my_pop = [idx,Population_class,APD_class,ER_class];

    % if delta < (M/100)*N
    for i = 1:cluster_size
        rr = find(idx==i);
    %       B = idx(idx(:,1) == i,:);
    %     APD = APD_class
    %     if size (rr,1)>1
    %     delta = 0;
        if (delta < 0.05*N)
            need = [];
            for j = 1:length(rr)
                APD_new = APD_class(rr(j)).c;
                [APD_new,index] = min(APD_new);
                need = [need;[rr(j),index,APD_new]];
            end
            [~,index_final] = min(need(:,end));
            need = need(index_final,:);
            sub_class = need(1);
            sub_class_entry = need(2);
            my_pop = [my_pop;Population_class(sub_class).c(sub_class_entry,:)];
    %          [~,ind] = min(B(:,end-1));
        %      [~,ind] = max(B(:,end));
        else 
             need = [];
            for j = 1:length(rr)
                ER_new = ER_class(rr(j)).c;
                [ER_new,index] = max(ER_new);
                need = [need;[rr(j),index,ER_new]];
%                 need = [need;[rr(j),ER_new]];
            end
            [~,index_final] = max(need(:,end));
            need = need(index_final,:);
            sub_class = need(1);
            sub_class_entry = need(2);
            my_pop = [my_pop;Population_class(sub_class).c(sub_class_entry,:)];
%              my_pop = [my_pop;Population_class(index_final).c];

        end


    %     else
    %         
    %         my_pop = [my_pop];
    %     end


    end
end
end

function [class,Empty_ref,fill_ref,APD_class,ER_class,Population_class] = reference(fitness,V,refV,ER,Population,theta)
[N,M] = size(fitness);
VN = size(V, 1);
Empty_ref = [];
empty_rows = [];
fill_ref = [];
% APD = [];

Zmin = min(fitness,[],1);	%?????
% Zmax = max(fitness,[],1);	%?????
%Z0 = min(Z0, Zmin)
%Translation
fitness = (fitness - repmat(Zmin, [size(fitness,1) 1]));

    
%???????
clear class;
ufitness = fitness./repmat(sqrt(sum(fitness.^2,2)), [1 M]);
% uV = V./repmat(sqrt(sum(V.^2,2)), [1 M]);
% cosine = ufitness*V'; %calculate the cosine values between each solution and each vector
cosine = ufitness*V'; %calculate the cosine values between each solution and each vector
acosine = acos(cosine);

[maxc maxcidx] = max(cosine, [], 2);
class = struct('c', cell(1,VN)); %classification
for i = 1:N
    class(maxcidx(i)).c = [class(maxcidx(i)).c, i];
end;
jj = 1;
Population_class = struct;
APD_class = struct;
ER_class = struct;
for k = 1:VN
    
%     if(sum(V(k,:)) >= 0.9999999 && sum(V(k,:)) <= 1.0005)
%         class(k).c = [1:N]; 
%     end
    if(isempty(class(k).c))
        Empty_ref = [Empty_ref;V(k,:)];
        empty_rows = [empty_rows;k];
    else
        sub = class(k).c;
        subFunctionValue = fitness(sub,:);
        Population_class(jj).c = Population(sub,:);
        subacosine = acosine(sub, k);
        subacosine = subacosine/refV(k);
        D1 = sqrt(sum(subFunctionValue.^2,2));
        APD_class(jj).c = D1.*(1 + (theta)*(subacosine));
        ER_class(jj).c = ER(sub,:);
        fill_ref = [fill_ref;V(k,:)];
        jj = jj+1;
    end
end
end