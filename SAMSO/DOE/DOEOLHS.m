function LHS = DOEOLHS(npoints, ndv, method, varargin)
%% from SRGTSToolbox
%最优拉丁超立方
%Function DOEOLHS generates optimized Latin hypercube samples. It uses
%either the enhanced stochastic evolutionary algorithm (ESEA) proposed by
%Jin et al. (2005) or the genetic algorithm (GA) proposed by Bates et al.
%(2004) to solve the optimization problem.
%Each row of the design represents one run (or sample). Design variables
%are normalized so that the hypercube points take values between 0 and 1.
%Thus, for example:
%
%     P = DOEOLHS(NPOINTS, NDV, METHOD): generates an NPOINTS-by-NDV
%     matrix. NPOINTS is the number of points and NDV is the number of
%     variables. METHOD is a string that defines the optimization
%     algorithm:
%         - 'ESEA' uses the enhanced stochastic evolutionary algorithm
%           proposed by Jin et al. (2005)
%         - 'GA' uses the genetic algorithm (GA) proposed by Bates et al.
%           (2004)
%
%     P = DOEOLHS(NPOINTS, NDV, 'ESEA', MAXITER, MAXSTALLITER):
%     creates the Latin hypercube design with the ESEA set to run for
%     MAXITER iterations with at most MAXSTALLITER iterations without
%     improvement.
%
%     P = srgtsDOEOLHS(NPOINTS, NDV, 'GA', MAXITER, MAXSTALLITER, POPSIZE):
%     creates the Latin hypercube design with the GA set to run for
%     MAXITER iterations with at most MAXSTALLITER iterations without
%     improvement using POPSIZE individuals.
%
%Example:
%     % create a 14x2 design.
%     NPOINTS = 14;
%     NDV     = 2;
%
%     P = srgtsDOEOLHS(NPOINTS, NDV, 'ESEA')
% 
%     P =
% 
%     0.3846    0.7692
%     0.8462    0.8462
%     0.9231    0.2308
%     0.6154    0.3846
%     0.2308    0.0769
%     0.5385    0.9231
%     1.0000    0.5385
%     0.6923    0.6154
%     0.3077    0.4615
%     0.0769    0.3077
%     0.1538    1.0000
%          0    0.6923
%     0.4615    0.1538
%     0.7692         0
%
% 
%OBSERVATION:
%In both ESEA and GA strategies, the obtained experimental design might
%change from one run to another due to the random nature of the heuristic
%optimization techniques used to solve the optimization problem.
% 
%For further reference, see:
%Jin R, Chen W and Sudjianto A, "An efficient algorithm for constructing
%optimal design of computer experiments," Journal of Statistical Planning
%and Inference, Vol. 134, pp 268?87, 2005.
%
%Bates SJ, Sienz J, and Toropov VV, "Formulation of the optimal Latin
%hypercube design of experiments using a permutation genetic algorithm,"
%45th AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics and Materials
%Conference, Palm Springs, CA, 19?2 April 2004. AIAA-2004-2011.

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

switch method
    case 'ESEA'
        
        if nargin == 3
            maxiter      = 50;
            maxstalliter = 5;
        else
            maxiter      = varargin{1};
            maxstalliter = varargin{2};
        end
        
        LHS =ESEAOLHSdesign(npoints, ndv, maxiter, maxstalliter);
        
    case 'GA'
        if nargin == 3
            maxiter      = 50;
            maxstalliter = 20;
            popsize   = 10*ndv;
        else
            maxiter      = varargin{1};
            maxstalliter = varargin{2};
            popsize      = varargin{3};
        end
        
        LHS = GAOLHSdesign(npoints, ndv, maxiter, maxstalliter, popsize);

end

return
