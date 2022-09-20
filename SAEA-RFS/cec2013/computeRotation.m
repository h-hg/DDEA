
function [Q R] = computeRotation(D)
    A = randn(D);

    Q = zeros(D);
    R = zeros(D);
    for j = 1:D
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end

end

