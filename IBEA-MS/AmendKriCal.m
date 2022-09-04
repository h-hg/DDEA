function [Offspring,MSE] = AmendKriCal(Pop, KModel, Global, OriMSE)
    %% parameter set
    OriDec = decs(Pop);
    SumMSE = sum(OriMSE,2);
    OriObj = objs(Pop);
    OriIND = find(SumMSE==0);
    %%% Maintaining 0 of MSE values if there exists the individuals from
    %%% the initial offline data
    IniDec = OriDec(OriIND, :);
    IniObj = OriObj(OriIND, :);
    IniMSE = OriMSE(OriIND, :);
    OriDec(OriIND,:) = [];
    Decs   = OriDec;  
    [N, ~] = size(Decs);
    Objs   = zeros(N, Global.M);
    MSE    = zeros(N, Global.M);
    %% Predict the Objs and MSEs of parents by the Kriging models
    if ~isempty(Decs)
        %%% (Dace tool is hard to predict only one individual in one turn)
        if N==1
            N    = 2;
            Decs = [Decs;Decs];
            Objs = zeros(N, Global.M);
            MSE  = zeros(N, Global.M);
            for j=1:Global.M
                [Objs(:,j), MSE(:,j)] = predictor(Decs, KModel(j).md);
            end
            Decs = Decs(1,:); Objs = Objs(1,:); MSE = MSE(1,:);
        else
            for j=1:Global.M
                [Objs(:,j), MSE(:,j)] = predictor(Decs, KModel(j).md);
            end
        end
        Decs = [Decs;IniDec];
        Objs = [Objs;IniObj];
        MSE  = [MSE; IniMSE];
        Offspring = PopStruct(Decs, Objs);
    else
        Decs = IniDec; Objs = IniObj; MSE = IniMSE;
        Offspring = PopStruct(Decs, Objs);
    end
end