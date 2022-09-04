%   DACE (Design and Analysis of Computer Experiments) Toolbox
%   Version 2.5, September 4, 2002
%   Copyright (c) 2002 by Hans Bruun Nielsen and IMM. 
% 
%   Model construction
%     dacefit - Constrained non-linear least-squares fit of a given
%               correlation model to the provided data set and
%               regression model.
% 
%   Model prediction
%     predictor - Model predictor with mean squared error estimate.
% 
%   Regression functions
%     regpoly0  - Zero order polynomial.
%     regpoly1  - First order polynomial.
%     regpoly2  - Second order polynomial.
% 
%   Correlation functions
%     corrcubic     - Local support, cubic polynomial
%     correxp       - Exponential.
%     correxpg      - General exponential.
%     corrgauss     - Gaussian.
%     corrlin       - Local support, linear.
%     corrspherical - Local support, spherical.
%     corrspline    - Local support, cubic spline.
% 
%   Experimental Design
%     gridsamp - Points in a regular grid.
%     lhsamp   - Latin hypercube distributed random numbers.
%
%   Auxiliary functions
%     dsmerge  - Merge data for multiple design sites.
%     
%   Data files
%     data1.mat - Example data S and Y 

