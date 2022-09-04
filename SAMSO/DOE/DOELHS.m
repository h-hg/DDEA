function LH =DOELHS(npoints, ndv, niterations)
%% from SRGTSToolbox
%Function DOELHS generates a Latin hypercube sample using random
%search. Function srgtsDOELHS is especially suitable for very large
%designs, when the native MATLAB function run out of memory. This may
%depend on the MATLAB version and on the configuration of the machine used.
%This function is most useful after you try to run 'lhsdesign' without
%success.Each row of the design represents one point (or sample). Design
%variables are normalized so that the hypercube points take values between
%0 and 1. Thus, for example:
%
%     P = DOELHS(NPOINTS, NDV, NITERATIONS): generates an
%     NPOINTS-by-NDV matrix. NPOINTS is the number of points and NDV is the
%     number of variables. It uses NITERATIONS iterations of the maximum
%     inter-distance algorithm.
%
%Example:
%     % create a 1024x10 design.
%     NPOINTS     = 1024;
%     NDV         = 10;
%     NITERATIONS = 5;
%
%     P = srgtsDOELHS(NPOINTS, NDV, NITERATIONS); % it may take a while

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
% I know, it is not a decent optimization, but it all I can do now...

J = -Inf;
LH = [];
for c1 = 1 : niterations

    LHtrial = [];
    for c2 = 1 : ndv
        LHtrial = [LHtrial, randperm(npoints).'];
    end

    if niterations > 1
        Jnew = DOEMinDistCriterion(LHtrial,npoints);
        if Jnew > J
            J  = Jnew; %JnewµÄ×÷ÓÃ£¿
            LH = LHtrial;
        end
    else
        LH = LHtrial;
    end

end

LH = ScaleVariable(LH, ...
                        [ones(1,ndv); npoints*ones(1,ndv)], ...
                        [zeros(1,ndv); ones(1,ndv)]);

return
