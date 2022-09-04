% % %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ------------------GRIEWANK FUNCTION
function [y] = griewank(xx)
%
% INPUT:
%
% xx = [x1, x2, ..., xd]
% % %
d =  length(xx);
sum = 0;
prod = 1;
for ii = 1:d
	xi = xx(ii);
	sum = sum + xi^2/4000;
	prod = prod * cos(xi/sqrt(ii));
end
y = sum - prod + 1;
end