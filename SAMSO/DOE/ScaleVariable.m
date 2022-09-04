function [xstar] = ScaleVariable(x, fromDomain, toDomain)
%Function srgtsScaleVariable scales a vector from a specific domain to
%another. Thus, for example:
%
%     xstar = srgtsScaleVariable(x, FROMDOMAIN, TODOMAIN): maps the vector
%     x from the FROMDOMAIN domain to the TODOMAIN domain.
%
%Example:
%     physicalspace = [-5  0;  % lower bound
%                      10 15]; % upper bound
%
%     normalizedspace = [0  0; % lower bound
%                        1  1]; % upper bound
%
%     % create points in the normalized space
%     P = [0.0  0.0
%          0.5  0.5
%          1.0  1.0];
%
%     % map P to the physical space
%     X = srgtsScaleVariable(P, normalizedspace, physicalspace)
%
%     X =
% 
%    -5.0000         0
%     2.5000    7.5000
%    10.0000   15.0000
% 
%     % create points in the physical space
%     X = [2.5  0.0
%          -5   7.5
%          10   15];
%
%     % map X to the normalized space
%     P = srgtsScaleVariable(X, physicalspace, normalizedspace)
% 
%     P =
% 
%     0.5000         0
%     0.0000    0.5000
%     1.0000    1.0000

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
[npoints, ndv] = size(x);

rangeFrom = fromDomain(2,:) - fromDomain(1,:);
rangeTo   = toDomain(2,:)   - toDomain(1,:);

a = rangeTo./rangeFrom;
b = toDomain(2,:) - a.*fromDomain(2,:);

xstar  = zeros(npoints, ndv);

for counter01 = 1 : npoints
    
    % mapping
    xstar(counter01,:) = a.*x(counter01,:) + b;
    
    % make sure that points in the lower bound will not suffer with MATLAB
    % precision
    mask = x(counter01,:) == fromDomain(1,:); % detects if point has any variable in the lower boud
    xstar(counter01,:) = mask.*toDomain(1,:) + ...
        (~mask).*xstar(counter01,:);
    
    % make sure that points in the upper bound will not suffer with MATLAB
    % precision
    mask = x(counter01,:) == fromDomain(2,:); % detects if point has any variable in the upper boud
    xstar(counter01,:) = mask.*toDomain(2,:) + ...
        (~mask).*xstar(counter01,:);
    
end

return
