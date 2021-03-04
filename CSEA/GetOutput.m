%Character with two types of solutions, 0 or 1
function Output = GetOutput(PopObj,RefPoint)
    N = size(PopObj,1);
    Output = true(N,1);
    for i = 1 : size(RefPoint,1)
        Output = Output & any(PopObj<=repmat(RefPoint(i,:),N,1),2);
    end
end