% full_crossover
function [U] = full_crossover(F, hisx)
    [NP Dim] = size(hisx);
    pos0 = hisx(1,:);
    a = randperm(Dim);
    for i=1:Dim 
        k = a(i);
        pos = repmat(pos0,NP,1);
        pos(:,k) = hisx(:,k);
        fit = F(pos);
        [bestfit index] = min(fit);
        pos0(k) = hisx(index,k); 
    end
    U = pos0;
end