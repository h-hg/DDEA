function [y] = ACKLEY(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ps,D]=size(x);
    y=sum(x.^2,2);
    y=20-20.*exp(-0.2.*sqrt(y./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
    y= y';
end