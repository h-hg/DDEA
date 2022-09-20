
function s = computeShiftVector(D, lb, ub)
    s = zeros(D, 1);
    hw = (ub - lb) / 2.0;
    middle = lb + hw;
    for i=1:D
        s(i) = middle + randn * hw;
        while((s(i) < lb) || (s(i) > ub))
            s(i) = middle + randn * hw;
        end
    end
end
