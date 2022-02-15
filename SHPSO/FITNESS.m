%%%*********************************************************************************************%%%
%% Implementation of Surrogate-assisted Hierarchical Particle Swarm Optimization (SHPSO)
%% H. Yu, Y. Tan, J. Zeng, C. Sun, Y. Jin, Surrogate-assisted hierarchical 
%% particle swarm optimization, Information Sciences, 454-455 (2018) 59-72.
%%%*********************************************************************************************%%%
%% This paper and this code should be referenced whenever they are used to 
%% generate results for the user's own research. 
%%%*********************************************************************************************%%%
%% This matlab code was written by Haibo Yu
%% Please refer with all questions, comments, bug reports, etc. to tyustyuhaibo@126.com
% 
%% Test functions for category: 'FITNESS'

function [y] = FITNESS(xx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % -------------------ACKLEY FUNCTION
% % function [y] = ackley(xx, a, b, c)
% % INPUTS:
% %
% % xx = [x1, x2, ..., xd]
% % a = constant (optional), with default value 20
% % b = constant (optional), with default value 0.2
% % c = constant (optional), with default value 2*pi
%
% d = length(xx);
% if (nargin < 4)
%     c = 2*pi;
% end
% if (nargin < 3)
%     b = 0.2;
% end
% if (nargin < 2)
%     a = 20;
% end
% sum1 = 0;
% sum2 = 0;
% for ii = 1:d
% 	xi = xx(ii);
% 	sum1 = sum1 + xi^2;
% 	sum2 = sum2 + cos(c*xi);
% end
% term1 = -a * exp(-b*sqrt(sum1/d));
% term2 = -exp(sum2/d);
% y = term1 + term2 + a + exp(1);
% end
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ------------------GRIEWANK FUNCTION
% % function [y] = griewank(xx)
% %
% % INPUT:
% %
% % xx = [x1, x2, ..., xd]
% % % %
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % -------------------------ROSENBROCK FUNCTION
% function [y] = rosen(xx)
% %
% % INPUT:
% %
% % xx = [x1, x2, ..., xd]
% %
% d = length(xx);
% sum = 0;
% for ii = 1:(d-1)
% 	xi = xx(ii);
% 	xnext = xx(ii+1);
% 	new = 100*(xnext-xi^2)^2 + (xi-1)^2;
% 	sum = sum + new;
% end
% y = sum;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % -------------------------Ellipsoid FUNCTION
% function [y] = ellipsoid(xx)
% % INPUT:
% %
% % xx = [x1, x2, ..., xd]
% %
% d = length(xx);
% sum = 0;
% for ii = 1:d
% 	xi = xx(ii);
% 	sum = sum + ii * xi^2;
% end
% y = sum;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
