% -------------------------Shifted RASTRIGIN FUNCTION
function [y] = S_RASTRIGIN(x)
    [ps,D]=size(x);
    Xmin=-5.12; Xmax=5.12;
    o = (Xmax-Xmin)/6;
    x=x-repmat(o,ps,1);
    y=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
    y= y';
end