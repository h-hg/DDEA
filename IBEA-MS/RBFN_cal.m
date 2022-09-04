function Offspring = RBFN_cal(Decs, obj_nets, Global)
    %%
    [N, ~] = size(Decs);
    Objs = zeros(N, Global.M);
    %% calculate for evr model && return the average predict
    kernal_name = obj_nets.name;
    eval(['Z = ', kernal_name, '(Decs, obj_nets.centers, obj_nets.sigma);']);
    Z    = [Z, ones(size(Z, 1), 1)];
    Objs = Z*obj_nets.weight+Objs;
    %%
    Offspring = PopStruct(Decs, Objs);
end