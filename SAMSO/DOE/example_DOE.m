%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is just an example of how to use the SURROGATES Toolbox
% capabilities.
%
% It illustrates the use of the following functions:
%     - srgtsDOEFullFactorial.m
%     - srgtsDOELHS.m
%     - srgtsDOEOLHS.m
%     - srgtsDOETPLHS.m
%     - srgtsDOESubSample.m
%
% All parameters used in this example just illustrates how modify the
% default ones. They do not reflect any specific study in parameter
% selection.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% two design variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndv     = 2;
npoints = 14;
nlevels = 3; % used in the full-factorial design

% full-factorial design (FFD)
disp(sprintf('full-factorial design (FFD)...'));
FFD = DOEFullFactorial(ndv, nlevels);

% Latin hypercube sample (LHS)
disp(sprintf('Latin hypercube sample (LHS)...'));
niterations    = 5;
LHS = DOELHS(npoints, ndv, niterations);

% enhanced Latin hypercube design (TPLHS) 
disp(sprintf('enhanced Latin hypercube design (TPLHS)...'));
TPLHS = DOETPLHS(npoints, ndv);

% optimal Latin hypercube design (ESEA) 
disp(sprintf('optimal Latin hypercube design (ESEA) ...'));
OLH_ESEA = srgtsDOEOLHS(npoints, ndv, 'ESEA');

% optimal Latin hypercube design (GA) 
disp(sprintf('optimal Latin hypercube design (GA) ...'));
OLH_GA = DOEOLHS(npoints, ndv, 'GA');

% selecting points from a user-defined desing (SP)
disp(sprintf('selecting points from a user-defined desing...'));
SSmaxmin  = DOESubSample(LHS, npoints/2, 'MaxMin', 1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf(1)
subplot(2,3,1);
plot(FFD(:,1), FFD(:,2), 'o'); title('FFD'); axis square;

subplot(2,3,2);
plot(LHS(:,1), LHS(:,2), 'o'); title('LHS'); axis square;

subplot(2,3,3);
plot(TPLHS(:,1), TPLHS(:,2), 'o'); title('TPLHS'); axis square;

subplot(2,3,4);
plot(OLH_ESEA(:,1), OLH_ESEA(:,2), 'o'); title('OLH-ESEA'); axis square;

subplot(2,3,5);
plot(OLH_GA(:,1), OLH_GA(:,2), 'o'); title('OLH-GA'); axis square;

subplot(2,3,6);
plot(SSmaxmin(:,1), SSmaxmin(:,2), 'o'); title('SSmaxmin'); axis square;
