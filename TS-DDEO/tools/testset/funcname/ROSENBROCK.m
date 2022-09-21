% % -------------------------ROSENBROCK FUNCTION
function [y] = ROSENBROCK(x)
    [ps,D]=size(x);
    y=sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);
    y= y';
end