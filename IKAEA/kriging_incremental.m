function model_new = kriging_incremental(model,add_x,add_y,sample_x,sample_y,lower_bound,upper_bound)
% parameters of old model
theta = model.theta;
L = model.L;
% number of incremental samples
n = size(sample_x,1);
k = size(add_x,1);
X = (sample_x - lower_bound)./(upper_bound - lower_bound);
dX = (add_x - lower_bound)./(upper_bound - lower_bound);
% calculate A matrix and B matrix
temp1 = sum(X.^2.*theta,2)*ones(1,k);
temp2 = sum(dX.^2.*theta,2)*ones(1,n);
A = exp(-max(temp1 + temp2'-2.*(X.*theta)*dX',0));
temp = sum(dX.^2.*theta,2)*ones(1,k);
B = exp(-max(temp + temp'-2.*(dX.*theta)*dX',0));
% update the L matrix
L21 = (L\A)';
% add a small value to the diagonal to avoid numerical problems
L22 = chol(B - L21*L21'+eye(k).*1E-10,'lower');
L_new = [L,zeros(n,k);L21,L22];
% calculate mu and sigma2
one = ones(n+k,1);
Y_new = [sample_y;add_y];
mu_new = (one'*(L_new'\(L_new\Y_new)))/(one'*(L_new'\(L_new\one)));
sigma2_new = ((Y_new-mu_new)'*(L_new'\(L_new\(Y_new-mu_new))))/(n+k);
% output the results of the DACE model
model_new.theta = theta;
model_new.mu = mu_new;
model_new.sigma2 = sigma2_new;
model_new.L = L_new;
model_new.lnL = -(-0.5*n*log(sigma2_new)-sum(log(abs(diag(L_new)))));
end












