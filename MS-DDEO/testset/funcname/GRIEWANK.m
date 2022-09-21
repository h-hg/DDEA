% % ------------------GRIEWANK FUNCTION
function [y] = GRIEWANK(x)
    [ps,D]=size(x);
    y=1;
    for i=1:D
        y=y.*cos(x(:,i)./sqrt(i));
    end
    y=sum(x.^2,2)./4000-y+1;
    y= y';
end