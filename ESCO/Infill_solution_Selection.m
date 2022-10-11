function new = Infill_solution_Selection(Data,S,l)

new = [];
pool = [];
Ns = 1;
K = 3;
best = min(Data(:,end));
X = Data(:,1:end-1);
x = S(:,1:end-1);  y = S(:,end);
Eu_dis = pdist2(x,X);
min_dis = min(Eu_dis,[],2);
logi1 = min_dis > l;
logi2 = y < best;
pool = [pool;S(logi1 & logi2,:)];

if ~isempty(pool)
    K = min(K,size(pool,1));
    [~,ind] = sort(pool(:,end));
    rng('shuffle');
    r = randperm(K,min(Ns,K));
    new = pool(ind(r),:);
elseif sum(logi1)
    pool2 = [pool;S(logi1,:)];
    [~,ind] = sort(pool2(:,end));
    rng('shuffle');
    r = randperm(size(pool2,1),1);
    new = pool2(ind(r),:);
else
    [~,ind] = sort(S(:,end));
    rng('shuffle');
    r = randperm(size(S,1),1);
    new = S(ind(r),:);
end

end

