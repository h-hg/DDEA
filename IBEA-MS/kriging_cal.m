function [Offspring,MSE] = kriging_cal(Decs, Mds, Global)
    
    [N, ~] = size(Decs);
    Objs   = zeros(N, Global.M);
    MSE    = zeros(N, Global.M);

    for j=1:Global.M
        [Objs(:,j), MSE(:,j)] = predictor(Decs, Mds(j).md);
    end
    Offspring = PopStruct(Decs, Objs);
end