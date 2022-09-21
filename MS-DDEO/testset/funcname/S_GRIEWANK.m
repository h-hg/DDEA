% % ------------------Shifted GRIEWANK FUNCTION
function [y] = S_GRIEWANK(x)
    [ps,D]=size(x);
    Xmin=-600; Xmax=600;
    o = (Xmax-Xmin)/6;
    x=x-repmat(o,ps,1);
    y=1;
    for i=1:D
        y=y.*cos(x(:,i)./sqrt(i));
    end
    y=sum(x.^2,2)./4000-y+1;
    y= y';
end