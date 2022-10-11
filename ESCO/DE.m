function x = DE(pop,y,bu,bd,F,CR,strategy)

%------Evaluate the best member after initialization----------------------
% ibest   = 1;                      % start with first population member
% val(1)  = feval(fname,pop(ibest,:),y);
% bestval = val(1);                 % best objective function value so far
% nfeval  = nfeval + 1;
% for i=2:NP                        % check the remaining members
%     val(i) = feval(fname,pop(i,:),y);
%     nfeval  = nfeval + 1;
%     if (val(i) < bestval)           % if member is better
%         ibest   = i;                 % save its location
%         bestval = val(i);
%     end
% end
NP0 = size(pop,1);
[~,ind] = sort(y);
pop = pop(ind(1:ceil(NP0/2)),:); y = y(ind(1:ceil(NP0/2)));
[NP,D] = size(pop); XVmin = bd; XVmax = bu;
[~,ibest] = min(y);
bestmemit = pop(ibest,:);         % best member of current iteration

%------DE-Minimization---------------------------------------------
%------popold is the population which has to compete. It is--------
%------static through one iteration. pop is the newly--------------
%------emerging population.----------------------------------------

pm1 = zeros(NP,D);              % initialize population matrix 1
pm2 = zeros(NP,D);              % initialize population matrix 2
pm3 = zeros(NP,D);              % initialize population matrix 3
pm4 = zeros(NP,D);              % initialize population matrix 4
pm5 = zeros(NP,D);              % initialize population matrix 5
bm  = zeros(NP,D);              % initialize bestmember  matrix
ui  = zeros(NP,D);              % intermediate population of perturbed vectors
mui = zeros(NP,D);              % mask for intermediate population
mpo = zeros(NP,D);              % mask for old population
rot = (0:1:NP-1);               % rotating index array (size NP)
rotd= (0:1:D-1);                % rotating index array (size D)
rt  = zeros(NP);                % another rotating index array
rtd = zeros(D);                 % rotating index array for exponential crossover
a1  = zeros(NP);                % index array
a2  = zeros(NP);                % index array
a3  = zeros(NP);                % index array
a4  = zeros(NP);                % index array
a5  = zeros(NP);                % index array
ind = zeros(4);

XVmin = repmat(XVmin, NP,1);
XVmax = repmat(XVmax, NP,1);

popold = pop;                   % save the old population

ind = randperm(4);              % index pointer array

a1  = randperm(NP);             % shuffle locations of vectors
rt = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
a2  = a1(rt+1);                 % rotate vector locations
rt = rem(rot+ind(2),NP);
a3  = a2(rt+1);
rt = rem(rot+ind(3),NP);
a4  = a3(rt+1);
rt = rem(rot+ind(4),NP);
a5  = a4(rt+1);

pm1 = popold(a1,:);             % shuffled population 1
pm2 = popold(a2,:);             % shuffled population 2
pm3 = popold(a3,:);             % shuffled population 3
pm4 = popold(a4,:);             % shuffled population 4
pm5 = popold(a5,:);             % shuffled population 5

for i=1:NP                      % population filled with the best member
    bm(i,:) = bestmemit;          % of the last iteration
end

mui = rand(NP,D) < CR;          % all random numbers < CR are 1, 0 otherwise

if (strategy > 5)
    st = strategy-5;		  % binomial crossover
else
    st = strategy;		  % exponential crossover
    mui=sort(mui');	          % transpose, collect 1's in each column
    for i=1:NP
        n=floor(rand*D);
        if n > 0
            rtd = rem(rotd+n,D);
            mui(:,i) = mui(rtd+1,i); %rotate column i by n
        end
    end
    mui = mui';			  % transpose back
end
mpo = mui < 0.5;                % inverse mask to mui

% strategy       1 --> DE/best/1/exp           6 --> DE/best/1/bin
%                2 --> DE/rand/1/exp           7 --> DE/rand/1/bin
%                3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
%                4 --> DE/best/2/exp           9 --> DE/best/2/bin
%                5 --> DE/rand/2/exp           else  DE/rand/2/bin

switch st
    case 1                      % DE/best/1
        ui = bm + F*(pm1 - pm2);    % differential variation
        ui = popold.*mpo + ui.*mui; % crossover
        
    case 2                      % DE/rand/1
        ui = pm3 + F*(pm1 - pm2);   % differential variation
        ui = popold.*mpo + ui.*mui; % crossover
        
    case 3                          % DE/rand-to-best/1
        ui = popold + F*(bm-popold) + F*(pm1 - pm2);
        ui = popold.*mpo + ui.*mui;     % crossover
        
    case 4                               % DE/best/2
        ui = bm + F*(pm1 - pm2 + pm3 - pm4); % differential variation
        ui = popold.*mpo + ui.*mui;          % crossover
        
    case 5                                % DE/rand/2
        ui = pm5 + F*(pm1 - pm2 + pm3 - pm4); % differential variation
        ui = popold.*mpo + ui.*mui;           % crossover
end

% correcting violations on the lower bounds of the variables
maskLB = ui > XVmin; % these are good to go
maskUB = ui < XVmax; % these are good to go
ui     = ui.*maskLB.*maskUB + XVmin.*(~maskLB) + XVmax.*(~maskUB);
x = ui;
end

