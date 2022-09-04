function model = kriging_theta_train(sample_x,sample_y,lower_bound,upper_bound,theta0,theta_lower,theta_upper)
[n,num_vari] = size(sample_x);
X = (sample_x - lower_bound)./(upper_bound - lower_bound);
Y = sample_y;
% optimize the theta values with in [10^a,10^b]
theta0 = log10(theta0);
theta_lower = log10(theta_lower);
theta_upper = log10(theta_upper);
options = optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',2*num_vari,'OptimalityTolerance',1E-20,'StepTolerance',1E-20,'Display','off');
[x_star,lnL] = fmincon(@(x)Simple_Concentrated_lnLikelihood(x,X,Y),theta0(1),[],[],[],[],theta_lower(1),theta_upper(1),[],options);
theta = 10.^(x_star*ones(1,num_vari));
one = ones(n,1);
% calculate the correlation matrix
temp1 = sum(X.^2.*theta,2)*one';
temp2 = X.*sqrt(theta);
R = exp(-(temp1 + temp1'-2.*(temp2*temp2'))) + eye(n).*(10+n)*eps;
% use the Cholesky factorization
L = chol(R,'lower');
% calculate mu and sigma
mu = (one'*(L'\(L\Y)))/(one'*(L'\(L\one)));
sigma2 = ((Y-mu)'*(L'\(L\(Y-mu))))/n;
% output the results of the DACE model
model.theta = theta;
model.mu = mu;
model.sigma2 = sigma2;
model.L = L;
model.lnL = lnL;

end




function  obj = Simple_Concentrated_lnLikelihood(x,X,Y)
num_vari = size(X,2);
theta = x*ones(1,num_vari);
theta = 10.^(theta);
% the concentrated ln-likelihood function
n = size(X,1);
one = ones(n,1);
% calculate the correlation matrix
temp1 = sum(X.^2.*theta,2)*one';
temp2 = X.*sqrt(theta);
R = exp(-(temp1 + temp1'-2.*(temp2*temp2'))) + eye(n).*(10+n)*eps;
% use the  Cholesky factorization
[L,p] = chol(R,'lower');
if p>0
    lnL = -1e8;
else
    mu = (one'*(L'\(L\Y)))/(one'*(L'\(L\one)));
    sigma2 = ((Y-mu)'*(L'\(L\(Y-mu))))/n;
    lnL = -0.5*n*log(sigma2)-sum(log(abs(diag(L))));
end
obj = -lnL;
end











