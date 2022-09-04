%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % -------------------ACKLEY FUNCTION
function [y] = ackley(xx, a, b, c)
% % INPUTS:
% %
% % xx = [x1, x2, ..., xd]
% % a = constant (optional), with default value 20
% % b = constant (optional), with default value 0.2
% % c = constant (optional), with default value 2*pi
% % %
d = length(xx);
if (nargin < 4)
    c = 2*pi;
end
if (nargin < 3)
    b = 0.2;
end
if (nargin < 2)
    a = 20;
end
sum1 = 0;
sum2 = 0;
for ii = 1:d
	xi = xx(ii);
	sum1 = sum1 + xi^2;
	sum2 = sum2 + cos(c*xi);
end
term1 = -a * exp(-b*sqrt(sum1/d));
term2 = -exp(sum2/d);
y = term1 + term2 + a + exp(1);
end