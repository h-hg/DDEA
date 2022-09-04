function X = GAOLHSdesign(npoints, ndv, maxiter, maxstalliter, popsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbfneval = popsize;
pop      = zeros(popsize, npoints*ndv);
scores   = zeros(popsize, 1);
for c1 = 1 : popsize
    ind          = encoding(createlhs(npoints, ndv), npoints, ndv);
    scores(c1,1) = myobjective(ind, npoints, ndv);
    pop(c1,:)  = ind;
end

phip = min(scores);
idx  = find(scores == phip);
P = pop(idx(1),:);

timetostop = 0; iter = 0; stalliter = 0;
while timetostop == 0

    iter = iter + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % selection
    [pop, scores] = OPTMGASelectionTournament(pop, scores, 2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % xover
    [pop, scores, nbfneval] = cycleXOver(pop, scores, popsize, npoints, ndv, nbfneval);
    [pop, scores, nbfneval] = inversionXOver(pop, scores, popsize, npoints, ndv, nbfneval);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mutation
    [pop, scores, nbfneval] = mutation(pop, scores, popsize, npoints, ndv, nbfneval);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time to stop?
    minscores = min(scores);
    if minscores < phip
        phip  = minscores;
        stalliter = 0;
        idxind = find(scores == minscores);
        P = pop(idxind(1),:);
    else
        stalliter = stalliter + 1;
    end

    timetostop = (iter == maxiter) || (stalliter == maxstalliter);

end

P = decoding(P, npoints, ndv);
if checkLHS(P, npoints, ndv)
    X = srgtsScaleVariable(P, ...
        [ones(1, ndv); npoints*ones(1, ndv)], ...
        [zeros(1, ndv); ones(1, ndv)]);
else
    X = [];
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friend functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = checkLHS(X, npoints, ndv)

v = 1;
idx = 1;
while ( (v == 1) && (idx <= ndv) )
    v = isequal(sort(X(:,idx)), [1:npoints]');
    idx = idx + 1;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = createlhs(npoints, ndv)

P = zeros(npoints,ndv);
for c1 = 1 : ndv
    P(:,c1) = randperm(npoints)';
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = encoding(X, npoints, ndv)

X = reshape(X,1,npoints*ndv);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = decoding(X, npoints, ndv)

X = reshape(X, npoints, ndv);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phip = myobjective(X, npoints, ndv)

phip = DOEPHIpCriterion(decoding(X, npoints, ndv));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newpop, newscores, nbfneval] = cycleXOver(pop, scores, popsize, npoints, ndv, nbfneval)

newpop    = pop;
newscores = scores;

children        = zeros(popsize, npoints*ndv);
childrenscores  = zeros(popsize,1);

p1idx = 1:2:popsize;
p2idx = 2:2:popsize;
for c1 = 1 : popsize/2

    p1 = decoding(pop(p1idx(c1), :), npoints, ndv);
    p2 = decoding(pop(p2idx(c1), :), npoints, ndv);

    ch1 = zeros(npoints, ndv);
    for c2 = 1 : ndv
        idx = 1; v = p1(idx,c2); idxrm = [2:npoints];
        flag = 0;
        while flag == 0
            ch1(idx,c2) = v;
            v = p2(idx,c2);
            if any(ch1(:,c2) == v)
                flag = 1;
            else
                idx = find(p1(:,c2) == v);
                idxrm(find(idxrm == idx)) = [];
            end

        end
        ch1(idxrm,c2) = p2(idxrm,c2);
    end

    ch2 = zeros(npoints, ndv);
    for c2 = 1 : ndv
        idx = 1; v = p2(idx,c2); idxrm = [2:npoints];
        flag = 0;
        while flag == 0
            ch2(idx,c2) = v;
            v = p1(idx,c2);
            if any(ch2(:,c2) == v)
                flag = 1;
            else
                idx = find(p2(:,c2) == v);
                idxrm(find(idxrm == idx)) = [];
            end

        end
        ch2(idxrm,c2) = p1(idxrm,c2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    children(p1idx(c1),:) = encoding(ch1, npoints, ndv);
    children(p2idx(c1),:) = encoding(ch2, npoints, ndv);
    if isequal(ch1, p1)
        childrenscores(p1idx(c1),:) = scores(p1idx(c1));
    elseif isequal(ch1, p2)
        childrenscores(p1idx(c1),:) = scores(p2idx(c1));
    else
        nbfneval = nbfneval + 1;
        childrenscores(p1idx(c1),:) = myobjective(encoding(ch1, npoints, ndv), npoints, ndv);
    end

    if childrenscores(p1idx(c1),:) < scores(p1idx(c1))
        newpop(p1idx(c1),:) = children(p1idx(c1),:);
        newscores(p1idx(c1),:) = childrenscores(p1idx(c1),:);
    end

    if isequal(ch2, p1)
        childrenscores(p2idx(c1),:) = scores(p1idx(c1));
    elseif isequal(ch2, p2)
        childrenscores(p2idx(c1),:) = scores(p2idx(c1));
    else
        nbfneval = nbfneval + 1;
        childrenscores(p2idx(c1),:) = myobjective(encoding(ch1, npoints, ndv), npoints, ndv);
    end

    if childrenscores(p2idx(c1),:) < scores(p2idx(c1))
        newpop(p2idx(c1),:) = children(p2idx(c1),:);
        newscores(p2idx(c1),:) = childrenscores(p2idx(c1),:);
    end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newpop, newscores, nbfneval] = inversionXOver(pop, scores, popsize, npoints, ndv, nbfneval)

newpop    = pop;
newscores = scores;

children        = zeros(popsize, npoints*ndv);
childrenscores  = zeros(popsize,1);

p1idx = 1:2:popsize;
p2idx = 2:2:popsize;
for c1 = 1 : popsize/2

    p1 = decoding(pop(p1idx(c1), :), npoints, ndv);
    p2 = decoding(pop(p2idx(c1), :), npoints, ndv);

    ch1 = zeros(npoints, ndv); ch2 = zeros(npoints, ndv);
    for c2 = 1 : ndv

        idx = sort(round(rand(2,1) * (npoints - 2)) + 1);

        if idx(1) == idx(2)
            ch1(:,c2) = p1(:,c2);
            ch2(:,c2) = p2(:,c2);
        else
            idxaux = [idx(1) + 1:idx(2)]; idxaux = sort(idxaux, 'descend');
            ch1(:,c2) = [p1(1:idx(1),c2)
                p1(idxaux,c2)
                p1(idx(2) + 1 : npoints,c2)];

            ch2(:,c2) = [p2(1:idx(1),c2)
                p2(idxaux,c2)
                p2(idx(2) + 1 : npoints,c2)];
        end
    end

    children(p1idx(c1),:) = encoding(ch1, npoints, ndv);
    children(p2idx(c1),:) = encoding(ch2, npoints, ndv);
    if isequal(ch1, p1)
        childrenscores(p1idx(c1),:) = scores(p1idx(c1));
    elseif isequal(ch1, p2)
        childrenscores(p1idx(c1),:) = scores(p2idx(c1));
    else
        nbfneval = nbfneval + 1;
        childrenscores(p1idx(c1),:) = myobjective(encoding(ch1, npoints, ndv), npoints, ndv);
    end

    if childrenscores(p1idx(c1),:) < scores(p1idx(c1))
        newpop(p1idx(c1),:) = children(p1idx(c1),:);
        newscores(p1idx(c1),:) = childrenscores(p1idx(c1),:);
    end

    if isequal(ch2, p1)
        childrenscores(p2idx(c1),:) = scores(p1idx(c1));
    elseif isequal(ch2, p2)
        childrenscores(p2idx(c1),:) = scores(p2idx(c1));
    else
        nbfneval = nbfneval + 1;
        childrenscores(p2idx(c1),:) = myobjective(encoding(ch1, npoints, ndv), npoints, ndv);
    end

    if childrenscores(p2idx(c1),:) < scores(p2idx(c1))
        newpop(p2idx(c1),:) = children(p2idx(c1),:);
        newscores(p2idx(c1),:) = childrenscores(p2idx(c1),:);
    end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newpop, newscores, nbfneval] = mutation(pop, scores, popsize, npoints, ndv, nbfneval)

newpop    = pop;
newscores = scores;

mutrate = popsize/10;
indidx  = round(rand(mutrate,1) * (popsize - 2)) + 1;
for c1 = 1 : mutrate
    newind = decoding(newpop(indidx(c1),:), npoints, ndv);
    for c2 = 1 : ndv
        varidx = round(rand(2,1) * (npoints - 2)) + 1;

        aux = newind(varidx(1),c2);
        newind(varidx(1),c2) = newind(varidx(2),c2);
        newind(varidx(2),c2) = aux;
    end
    newpop(indidx(c1),:)    = encoding(newind, npoints, ndv);
    newscores(indidx(c1),:) = myobjective(newind, npoints, ndv);
    nbfneval = nbfneval + 1;
end

return
