function Z = the_gaussian(PopDec, centers, sigma)
    r = pdist2(PopDec, centers);
    Z = exp(-(r.^2)./(2*sigma.^2));
end