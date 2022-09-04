function Rnets = construct_Rnets(pop, Global)

    PopObj     = objs(pop);
    center_num = ceil(sqrt(Global.D*11-1));   
    RName      = 'the_gaussian';
    Rnets.centers= get_center(pop, Global, center_num);          %Calculate the center points
    Rnets.sigma  = max(pdist(Rnets.centers))*5;
    Rnets.name   = RName;
    Z = get_Z(pop, RName, Rnets, Global);                        
    Rnets.weight = inv(Z'*Z)*Z'*PopObj;
    
end