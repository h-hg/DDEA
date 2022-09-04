function NbCoeff = srgtsPRSNumberOfCoefficients(NDV,PRSdegree)
%Function srgtsPRSNumberOfCoefficients calculates the number of
%coefficients in a PRS. Thus, for example:
%
%     NBCOEFF = srgtsPRSNumberOfCoefficients(NDV,PRSDEGREE): calculates the
%     number of coefficients in a PRS of degree PRSDEGREE for a model with
%     NDV variables.
% 
%                 (NDV + PRSDEGREE)!
%     NbCoeff = ----------------------
%               (NDV!) x (PRSDEGREE!)
%
%Example:
%     % model details
%     NDV        = 6;
%     PRSdegree  = 2;
%
%     % number of coefficients of the PRS
%     NbCoeff = srgtsPRSNumberOfCoefficients(NDV,PRSdegree)
% 
%     NbCoeff =
% 
%     28

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
% define the number of terms
NbCoeff = factorial(NDV + PRSdegree)/(factorial(NDV)*factorial(PRSdegree));

return
