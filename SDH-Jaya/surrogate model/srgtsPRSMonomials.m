function monomials = srgtsPRSMonomials(inp1, inp2)
%Function srgtsPRSMonomials outputs the vector of monomial used in the
%polynomial. Numbers in the string represent the index of individual
%variables. The constant term is represented by a zero. Variables in a
%monomial are separated by dots.
%
%     MONOMIALS = srgtsPRSMonomials(NDV, PRS_Degree): creates MONOMIALS
%     for the given number of variables NDV and PRS_Degree.
%
%     MONOMIALS = srgtsPRSMonomials(options, surrogate): creates MONOMIALS
%     vector to be used in conjunction with PRS_Beta in surrogate structure
%     to determine the PRS model fitted to data. Each beta coefficient
%     corresponds with a monomial. It also works for structures with
%     multiple surrogates; each surrogate's monomials are in a separate
%     column in the cell array.
%
%Example 1:
%     % basic information about the problem
%     ndv        = 3;
%     PRS_Degree = 2;
%
%     monomials = srgtsPRSMonomials(ndv, PRS_Degree)
%
%     monomials =
%
%     '0'       <- this is the constant term (intercept)
%     '1'       <- x1 term
%     '2'       <- x2 term
%     '3'       <- x3 term
%     '1.1'     <- x1*x1 term
%     '1.2'     <- x1*x2 term
%     '1.3'     <- x1*x3 term
%     '2.2'     <- x2*x1 term
%     '2.3'     <- x2*x3 term
%     '3.3'     <- x3*x3 term
% 
% 
%Example 2:
%     % basic information about the problem
%     myFN = @cos;  % this could be any user-defined function
%     designspace = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     % create DOE
%     npoints = 5;
%     X = linspace(designspace(1), designspace(2), npoints)';
%
%     % evaluate analysis function at X points
%     Y = feval(myFN, X);
%
%     % fit surrogate models
%     options   = srgtsPRSSetOptions(X, Y);
%     surrogate = srgtsPRSFit(options);
%
%     monomials = srgtsPRSMonomials(options, surrogate)

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
if isstruct(inp1)
    NbVariables    = inp2.NbVariables;
    PRS_Degree     = inp1.PRS_Degree;
    PRS_RemovedIdx = inp2.PRS_RemovedIdx;
else
    NbVariables    = inp1;
    PRS_Degree     = inp2;
    PRS_RemovedIdx = [];
end

m = {'1'};
for c1 = 1 : NbVariables
    m{c1,1} = num2str(c1);
end
monomials = ['0'; m];

for c1 = 1:PRS_Degree-1
    m2 = {};
    ci = 1;
    for c2 = 1:NbVariables
        while ~strncmp(m{ci},num2str(c2),floor(log10(c2))+1)
            ci = ci+1;
        end
        for c3 = ci:length(m)
            m2 = [m2;strcat(num2str(c2),'.',m{c3})];
        end
    end
    monomials = [monomials;m2];
    m = m2;
end

monomials(PRS_RemovedIdx) = [];

return