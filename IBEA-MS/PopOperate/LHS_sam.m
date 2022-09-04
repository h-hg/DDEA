function IniDec = LHS_sam(Global)
    N  = Global.D*11-1;
    IniDec = zeros(N, Global.D);
    for j=1:Global.D
        linepart = linspace(Global.lower(j), Global.upper(j), N+1)';
        randval  = rand(N,1)*(linepart(2)-linepart(1));
        IniDec(:,j) = linepart(1:end-1)+randval;
    end

    for j=1:Global.D
        ran_index    = randperm(size(IniDec, 1));
        IniDec(:, j) = IniDec(ran_index, j);
    end
end