function Offspring = AmendRBFCal(Pop, Rnets, Global, OriMSE)
    OriDec = decs(Pop);
    OriMSE = sum(OriMSE,2);
    OriObj = objs(Pop);
    OriIND = find(OriMSE==0);
    IniDec = OriDec(OriIND, :);
    IniObj = OriObj(OriIND, :);
    OriDec(OriIND,:) = [];
    Decs   = OriDec;  
    [N, ~] = size(Decs);
    Objs = zeros(N, Global.M);
    %% calculate objs by RBFNs
    kernal_name = Rnets.name;
    eval(['Z = ', kernal_name, '(Decs, Rnets.centers, Rnets.sigma);']);
    Z    = [Z, ones(size(Z, 1), 1)];
    Objs = Z*Rnets.weight+Objs;

    Decs = [Decs;IniDec];
    Objs = [Objs;IniObj];
    Offspring = PopStruct(Decs, Objs);
end