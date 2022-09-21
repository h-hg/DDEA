% -------------------------Ellipsoid FUNCTION
function [y] = Ellipsoid(x)
    [ps,d] = size(x);
    for i=1:d
        y=sum(i*x.^2,2);
    end
    y= y';
end

