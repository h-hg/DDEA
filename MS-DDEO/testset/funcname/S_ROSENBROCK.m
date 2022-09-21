% % -------------------------Shifted ROSENBROCK FUNCTION
function [y] = S_ROSENBROCK(x)
    [ps,D]=size(x);
    Xmin=-2.048; Xmax=2.048;
    o = (Xmax-Xmin)/6;
    x=x-repmat(o,ps,1);
    y=sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);
    y= y';
end