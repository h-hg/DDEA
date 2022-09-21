function U = Limited_boundary(X, lb, ub)
global Dim
X = X(:, 1:Dim);
[NP,~] = size(X);

for i = 1:NP
    X(i,:) = max(X(i,:), lb(1,:));
    X(i,:) = min(X(i,:), ub(1,:));
end
U = X;