function Population = Fitness (PopDec,Global)

    switch Global.problem
        case 'DTLZ1'
            g   = 100*(Global.D-Global.M+1+sum((PopDec(:,Global.M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,Global.M:end)-0.5)),2));
            Obj = 0.5*repmat(1+g,1,Global.M).*fliplr(cumprod([ones(size(PopDec,1),1),PopDec(:,1:Global.M-1)],2)).*[ones(size(PopDec,1),1),1-PopDec(:,Global.M-1:-1:1)];       
    end
    Population = PopStruct(PopDec, Obj);
end