function centers = get_center(pop, Global, center_num)
    PopDec  = decs(pop);
    [N, D]  = size(PopDec);
    [~,centers] = kmeans(PopDec,center_num);
end