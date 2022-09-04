function FFD =DOEFullFactorial(varargin)
%% from SRGTSToolbox
%Function DOEFullFactorial generates a mixed-level full-factorial   
%design (FFD).
% 
% THIS FUNCTION USES THE STATISTICS TOOLBOX!
% 
%Each row of the design represents one point (or sample). Design variables
%are normalized so that the hypercube points take values between 0 and 1.
%Thus, for example:
%
%     P = DOEFullFactorial(NDV,NLEVELS): generates an NPOINTS-by-NDV
%     matrix. NPOINTS is the number of points and NDV is the number of
%     variables.
%     NPOINTS = NLEVELS^NDV.
%
%     P = DOEFullFactorial(LEVELS): generates an NPOINTS-by-NDV
%     matrix. NPOINTS is the number of points and NDV is the number of
%     variables. LEVELS is an 1xNDV vector that specifies the number levels
%     for each design variable.
%
%Example 1:
%     nlevels = 4;
%     ndv     = 2;
%
%     P = DOEFullFactorial(ndv, nlevels)
% 
%     P =
% 
%     0.0000    0.0000
%     0.3333    0.0000
%     0.6667    0.0000
%     1.0000    0.0000
%     0.0000    0.3333
%     0.3333    0.3333
%     0.6667    0.3333
%     1.0000    0.3333
%     0.0000    0.6667
%     0.3333    0.6667
%     0.6667    0.6667
%     1.0000    0.6667
%     0.0000    1.0000
%     0.3333    1.0000
%     0.6667    1.0000
%     1.0000    1.0000
%
%     This generates a 16 point design with 4 levels for both of the two
%     design variables, i.e., it creates a 16x2 design.
%
%Example 2:
%     levels = [2 4 3];
%
%     P =DOEFullFactorial(levels)
% 
%     P =
% 
%          0    0.0000         0
%     1.0000    0.0000         0
%          0    0.3333         0
%     1.0000    0.3333         0
%          0    0.6667         0
%     1.0000    0.6667         0
%          0    1.0000         0
%     1.0000    1.0000         0
%          0    0.0000    0.5000
%     1.0000    0.0000    0.5000
%          0    0.3333    0.5000
%     1.0000    0.3333    0.5000
%          0    0.6667    0.5000
%     1.0000    0.6667    0.5000
%          0    1.0000    0.5000
%     1.0000    1.0000    0.5000
%          0    0.0000    1.0000
%     1.0000    0.0000    1.0000
%          0    0.3333    1.0000
%     1.0000    0.3333    1.0000
%          0    0.6667    1.0000
%     1.0000    0.6667    1.0000
%          0    1.0000    1.0000
%     1.0000    1.0000    1.0000
%
% This generates a 24 point design with 2 levels in the first design
% variable (0 and 1), 4 in the second one (0, 0.333, 0.667, and 1), and 3
% in the third one (0, 0.5, and 1). It creates a 24x3 design.

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

% check inputs
switch nargin
    case 1 % levels
        levels = varargin{1};
        ndv    = length(levels);

    case 2 % ndv, nlevels
        % NbPoints = P^ndv;
        ndv     = varargin{1};
        nlevels = varargin{2};

        levels = nlevels*ones(1,ndv);

end

FFD = ScaleVariable(fullfact( levels ), ...
                         [ones(1,ndv); levels], ...
                         [zeros(1,ndv); ones(1,ndv)]);

return
