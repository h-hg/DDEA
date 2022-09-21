function [lb, ub] = get_boundary(X)
global Dim

[NP, ~] = size(X);
lb = zeros(1, Dim);
ub = zeros(1, Dim);
for i = 1:Dim
    ub(1, i) = max(X(:, i));
    lb(1, i) = min(X(:, i));
end
% lb = lb - 0.0001;
% ub = ub + 0.0001;
