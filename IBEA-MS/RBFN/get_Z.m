function Z = get_Z(pop, RName, obj_nets, Global)
    kernal_name = RName;
    PopDec = decs(pop);
    eval(['Z = ', kernal_name, '(PopDec, obj_nets.centers, obj_nets.sigma);']);
    Z = [Z, ones(size(Z, 1), 1)];
end