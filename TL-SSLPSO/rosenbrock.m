% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% -------------------------ROSENBROCK FUNCTION
function [y] = rosenbrock(xx)
%
% INPUT:
%
% xx = [x1, x2, ..., xd]
% %
d = length(xx);
sum = 0;
for ii = 1:(d-1)
	xi = xx(ii);
	xnext = xx(ii+1);
	new = 100*(xnext-xi^2)^2 + (xi-1)^2;
	sum = sum + new;
end
y = sum;
end