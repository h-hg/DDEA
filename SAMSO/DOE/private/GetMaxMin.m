function [SubDoE, SubDoEIdx, RemDoE, RemDoEIdx] =GetMaxMin(DoE,npointstotal,npoints, ndv, niterations)

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

SubDoE    = [];
SubDoEIdx = [];
RemDoE    = [];
RemDoEIdx = [];
    
ncandidates = factorial(npointstotal)/(factorial(npoints)*factorial(npointstotal - npoints));
niterations = min(niterations, ncandidates);

% I know, it is not a decent optimization, but it all I can do now...
J = -Inf;
for c1 = 1 : niterations
    
    aux = randperm(npointstotal);
    
    idx   = aux(1 : npoints);
    idxrp = aux(npoints + 1 : end);
    
    X = DoE( idx , : );
    Jnew = DOEMinDistCriterion(X,npoints); % proper even for huge designs
    
    if Jnew > J
        J = Jnew;
        
        SubDoEIdx = idx;
        RemDoEIdx     = idxrp;

        SubDoE      = DoE(SubDoEIdx, : );
        RemDoE = DoE(RemDoEIdx, : );
    end
    
end

return
