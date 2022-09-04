% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % -------------------------Ellipsoid FUNCTION
function [y] = ellipsoid(xx)
% INPUT:
%
% xx = [x1, x2, ..., xd]
% % % 
d = length(xx);
sum = 0;
for ii = 1:d
	xi = xx(ii);
	sum = sum + ii * xi^2;
end
y = sum;
end