function [Selection,APD,my_ref,Empty_ref] = F_select(FunctionValue, V, theta, refV)

my_ref = [];
Empty_ref=[];
APD = [];
[N,M] = size(FunctionValue);
VN = size(V, 1);

Zmin = min(FunctionValue,[],1);
Zmax = max(FunctionValue,[],1);

FunctionValue = (FunctionValue - repmat(Zmin, [size(FunctionValue,1) 1]));

clear class;
for i=1:N
    uFunctionValue(i,:) = FunctionValue(i,:)./norm(FunctionValue(i,:));
end
%uFunctionValue = FunctionValue./repmat(sqrt(sum(FunctionValue.^2,2)), [1 M]);
cosine = uFunctionValue*V'; %calculate the cosine values between each solution and each vector
acosine = acos(cosine);

[maxc maxcidx] = max(cosine, [], 2);
class = struct('c', cell(1,VN)); %classification
for i = 1:N
    class(maxcidx(i)).c = [class(maxcidx(i)).c, i];
end;

Selection = []; FitnessValue = [];
EmptySet = []; NonEmptySet = [];


for k = 1:VN
%      if(sum(V(k,:)) >= 0.9999999 && sum(V(k,:)) <= 1.0005)
%         class(k).c = [1:N]; 
%      end
    
    if(~isempty(class(k).c))
        sub = class(k).c;
        subFunctionValue = FunctionValue(sub,:);
        subacosine = acosine(sub, k);
        D1 = sqrt(sum(subFunctionValue.^2,2));
        subacosine = subacosine/refV(k);
      
            if(sum(V(k,:)) >= 0.9999999 && sum(V(k,:)) <= 1.0005)
            %the extreme solution on each axis must be selected
                D = D1.*(1 + 1e200*(theta)*(subacosine));
            else
                D = D1.*(1 + (theta)*(subacosine));
            end;
        
        [mind mindidx] = min(D);
        
      
         Selection = [Selection; sub(mindidx)];
         APD = [APD;mind];
         my_ref = [my_ref;V(k,:)];
    else
        EmptySet = [EmptySet; k];
        Empty_ref = [Empty_ref;V(k,:)];
    end;
end;

Class = class(NonEmptySet);
end

