function Population = PopStruct(PopDec, PopObj)
    if isempty(PopDec)
        Population = struct;
    else
        for i=1:size(PopDec,1)
            Population(i).dec = PopDec(i,:);
            Population(i).obj = PopObj(i,:);
        end
    end
end