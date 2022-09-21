% -------------------------Shifted ELLIPSOID FUNCTION
function [y] = S_ELLIPSOID(x)
    [ps,d] = size(x);
    Xmin=-5.12; Xmax=5.12;
    o = (Xmax-Xmin)/6;
    x=x-repmat(o,ps,1);
    for i=1:d
        y=sum(i*x.^2,2);
    end
    y= y';
end

