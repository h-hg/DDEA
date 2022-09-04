function f = Ellipsoid(x)
% the Ellipsoid function
% xi = [-5.12,5.12]
f = sum((1:size(x,2)).*x.^2,2);
end