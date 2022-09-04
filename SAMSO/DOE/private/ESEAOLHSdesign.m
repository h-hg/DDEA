function X = ESEAOLHSdesign(npoints, ndv, maxiter, maxstalliter)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial LHD
P = createlhs(npoints, ndv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize ESEA parameters
p = 50; t = 1;
[phip D] = myObjective(P, npoints, p, t); nbfneval = 1;
Pbest    = P;
phipBest = phip;
Th0      = 0.005*phip;
Th       = Th0;
tol      = 0.025;

if(npoints < 23)
    M = ceil(20/(npoints - 1));
else
    M = round(2*npoints/50);
end
J = min(ceil((npoints*npoints - npoints) / 10.0 ), 50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outer loop itself
timetostop = 0; iter = 0; stalliter = 0;
while timetostop == 0
    iter = iter + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % before inner loop
    phipOldBest = phipBest;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % inner loop
    [P, D, phip, Pbest, phipBest, acptRatio, impRatio, nbfneval] = innerLoop(P, D, phip, Pbest, phipBest, p, t, npoints, ndv, M, J, Th, nbfneval);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % after inner loop
    if ((phipOldBest - phipBest)/phipBest > tol)
        flagImp = 1;
        stalliter = 0;
    else
        stalliter = stalliter + 1;
        flagImp = 0;
    end

    Th = updateTh(Th, flagImp, acptRatio, impRatio);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time to stop?
    timetostop = isItTimeToStop(iter, stalliter, maxiter, maxstalliter);

end

if checkLHS(Pbest, npoints, ndv)
    X = srgtsScaleVariable(Pbest, ...
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
function [P, D, phip, Pbest, phipBest, acptRatio, impRatio, nbfneval] = innerLoop(P, D, phip, Pbest, phipBest, p, t, npoints, ndv, M, J, Th, nbfneval)

% initialize parameters
ctrl  = 0;
nacpt = 0;
nimp  = 0;

% loop itself
while ctrl < M
    ctrl = ctrl + 1;

    % randomly pick J distinct element-exchanges within column
    % (ctrl mod ndv). Then, choose the best design Ptry from J designs
    % induced by exchanges
    col  = mod(ctrl, ndv) + 1;
    Ptry = P;
    idx = randperm(npoints); idx(3:end) = [];
    aux = Ptry(idx(1), col);
    Ptry(idx(1), col) = Ptry(idx(2), col);
    Ptry(idx(2), col) = aux;
    [phipTry Dtry] = myObjective(Ptry, npoints, p, t, D, col, idx);
    for c1 = 2 : J
        Paux = P;
        idx = randperm(npoints); idx(3:end) = [];
        aux = Paux(idx(1), col);
        Paux(idx(1), col) = Paux(idx(2), col);
        Paux(idx(2), col) = aux;
        [phipAux Daux] = myObjective(Paux, npoints, p, t, D, col, idx);
        if phipAux < phipTry
            Ptry    = Paux;
            phipTry = phipAux;
            Dtry    = Daux;
        end
    end
    nbfneval = nbfneval + J;

    % check acceptance and improvement
    if ((phipTry - phip) <= Th*rand)
        P     = Ptry;
        phip  = phipTry;
        D     = Dtry;
        nacpt = nacpt + 1;
        if phip < phipBest
            Pbest    = P;
            phipBest = phip;
            nimp = nimp + 1;
        end
    end

end

acptRatio = nacpt/M;
impRatio  = nimp/M;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Th = updateTh(Th, flagImp, acptRatio, impRatio)

aux = 1.0/0.7;
if(flagImp)
    alpha = 1.25;
    if ( (acptRatio > 0.1) && (impRatio < acptRatio) )
        alpha = 0.8;
    elseif ( (acptRatio > 0.1) && (impRatio == acptRatio) )
        alpha = 1.0;
    end

else
    if ( acptRatio < 0.1 )
        aux = 1.0/0.7;
    elseif ( acptRatio > 0.8 )
        aux = 0.9;
    end
    alpha = aux;
end

Th = alpha*Th;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = isItTimeToStop(iter, stalliter, maxiter, maxstalliter)

v = (iter == maxiter) || (stalliter == maxstalliter);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phip D] = myObjective(P, npoints, p, t, varargin)

if nargin == 4
    D = zeros(npoints,npoints);
    phip = 0;
    for c1 = 1 : npoints
        x01 = P(c1,:);
        for c2 = (c1 + 1) : npoints
            x02  = P(c2,:);
            d    = norm((x01 - x02), t);
            phip = phip + d^(-p);
            D(c1,c2) = d;
            D(c2,c1) = d;
        end
    end
else
    Dold = varargin{1};
    col  = varargin{2};
    idx  = varargin{3};
    elem01 = idx(1);
    elem02 = idx(2);

    D = updateDmatrix(Dold, P, npoints, col, elem01, elem02, t);

    phip = 0;
    for c1 = 1 : npoints
        for c2 = (c1 + 1) : npoints
            d = D(c1,c2);
            phip = phip + d^(-p);
        end

    end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = updateDmatrix(D, P, npoints, col, elem01, elem02, t)

for c1 = 1 : npoints
    if ( (c1 ~= elem01) && (c1 ~= elem02) )
        di1jnew = norm(P(c1,:) - P( elem01, : ),t);
        di2jnew = norm(P(c1,:) - P( elem02, : ),t);

        D(elem01,c1) = di1jnew;
        D(c1,elem01) = di1jnew;

        D(elem02,c1) = di2jnew;
        D(c1,elem02) = di2jnew;
    end
end

return
