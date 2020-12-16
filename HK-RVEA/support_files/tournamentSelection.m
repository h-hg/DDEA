function winners = tournamentSelection(tournSize, fitness)

n = length(fitness);
winners = [];

% Make sure that the population size is divisible by the tournament size
if(mod(n,tournSize) > 0)
    error('Population size has to be divisible by the tournament size\n');
end

% Repeat the process "tournSize" times
for i = 1:tournSize
    
    %Create a random set of competitors
    shuffleOrder = randperm(n);
    competitors = reshape(shuffleOrder, tournSize, n/tournSize)';
    
    %The winner is the competitor with best fitness
    [winFit, winID] = max(fitness(competitors),[],2);
    idMap = (0:tournSize-1)*n/tournSize;
    idMap1 = idMap(winID) + (1:size(competitors,1));
    winners = [winners; competitors(idMap1)];
end