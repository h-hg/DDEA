function [u, s] = kriging_predictor(test_x,model,sample_x,sample_y,lower_bound,upper_bound)
% parameters of Kriging model
theta = model.theta;
mu = model.mu;
sigma2 = model.sigma2;
L = model.L;
% normalize data
if ~isempty(lower_bound)
    X = (sample_x - lower_bound)./(upper_bound - lower_bound);
    x = (test_x - lower_bound)./(upper_bound- lower_bound);
else
    X = sample_x;
    x = test_x;
end
% initialize the prediction and variance
one = ones(size(sample_x,1),1);
% point-wise calculation
temp1 = sum(x.^2.*theta,2)*ones(1,size(X,1));
temp2 = sum(X.^2.*theta,2)*ones(1,size(x,1));
R = exp(-(temp1 + temp2'-2.*(x.*theta)*X'))';
u = mu + R' *(L'\(L\(sample_y - mu)));
mse = sigma2*(1 + (1-one'*(L'\(L\R)))'.^2/(one'*(L'\(L\one))) - sum((L\R).^2,1)');
s = sqrt(max(mse,0));


