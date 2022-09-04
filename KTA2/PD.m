function score = PD(Population,~)
% <max>
% Pure diversity
%------------------------------- Reference --------------------------------
% H. Wang, Y. Jin, and X. Yao, Diversity assessment in many-objective
% optimization, IEEE Transactions on Cybernetics, 2017, 47(6): 1510-1522.

      PopObj = Population;   
    if isempty(PopObj)
        score = nan;
    else
        C = false(length(Population));
        C(logical(eye(size(C)))) = true;
        D = pdist2( PopObj, PopObj,'minkowski',0.1);
        D(logical(eye(size(D)))) = inf;
        score = 0;
        for k = 1 : length(Population)-1
            while true
                [d,J] = min(D,[],2);
                [~,i] = max(d);
                if D(J(i),i) ~= -inf
                    D(J(i),i) = inf;
                end
                if D(i,J(i)) ~= -inf
                    D(i,J(i)) = inf;
                end
                P = any(C(i,:),1);
                while ~P(J(i))
                    newP = any(C(P,:),1);
                    if P == newP
                        break;
                    else
                        P = newP;
                    end
                end
                if ~P(J(i))
                    break;
                end
            end
            C(i,J(i)) = true;
            C(J(i),i) = true;
            D(i,:)    = -inf;
            score     = score + d(i);
        end
    end
end