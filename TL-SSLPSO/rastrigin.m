% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % --------------------------RASTRIGIN FUNCTION
function [y] = rastrigin(xx)
%
% INPUT:
%
% xx = [x1, x2, ..., xd]
% % %
d = length(xx);
sum = 0;
for ii = 1:d
	xi = xx(ii);
	sum = sum + (xi^2 - 10*cos(2*pi*xi));
end
y = 10*d + sum;
end