% % -------------------------Shifted ROSENBROCK FUNCTION
function [y] = S_ACKLEY(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ps,D]=size(x);
    Xmin=-32.768; Xmax=32.768;
    o = (Xmax-Xmin)/6;
    x=x-repmat(o,ps,1);
    y=sum(x.^2,2);
    y=20-20.*exp(-0.2.*sqrt(y./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
    y= y';
end