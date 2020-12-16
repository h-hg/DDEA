function [MatingPool] = F_mating(Population,pop_size)

    if size(Population,1) < pop_size
        N = pop_size;
        D = size(Population,2);
    else
        [N,D] = size(Population);
    end
    MatingPool = zeros(N,D);
    
    for i=1:pop_size
        RandList = randperm(size(Population,1));
        MatingPool(i,:) = Population(RandList(1), :);
    end
    if(mod(N,2) == 1)
        MatingPool = [MatingPool; MatingPool(1,:)];
    end;
end
