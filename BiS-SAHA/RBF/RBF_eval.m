function Yest=RBF_eval(X,S,lambda,gamma,flag)
% function Yest=RBF_eval(X,S,lambda,gamma,rho,flag)   % For Gaussian basis function

%--------------------------------------------------------------------------
%Copyright (c) 2012 by Juliane Mueller
%
% This file is part of the surrogate model module toolbox.
%
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%Tampere University of Technology, Finland
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%input:
%X are points where function values should be calculated
%S are points where the function values are known
%lambda parameter vector
%gamma contains the parameters of the optional polynomial tail
%flag is a string indicating which RBF to be used
%output: 
%the estimated function values at the points in X
%--------------------------------------------------------------------------

[mX,nX]=size(X); %dimensions of the points where function value should be predicted
[mS,nS]=size(S); %dimesnions of sample sites taken so far 
if nX~=nS %check that both matrices are of the right shape
    X=X';
    [mX,nX]=size(X);
end

% R=zeros(mX,mS); %compute pairwise distances of points in X and S
% for ii = 1:mX
%     for jj = 1:mS
%         R(ii,jj)=norm(X(ii,:)-S(jj,:));
%     end
% end

% R = single(R);
% X = single(X);
% S = single(S);
% for ii = 1:mX
%     XX = repmat(X(ii,:), mS, 1);
%     difference = XX - S;
%     square = difference.^2;
%     summation = sum(square, 2);
%     distance = realsqrt(summation);
%     R(ii, :) = distance';
% end

R = sqrt(X.^2*ones(size(S'))+ones(size(X))*(S').^2-2*X*S');



% for ii = 1:mX
%     XX = repmat(Xg(ii,:), mS, 1);
%     difference = XX - Sg;
%     square = difference.^2;
%     summation = sum(square, 2);
%     distance = realsqrt(summation);
%     R(ii, :) = distance';
% end
% R = gather(R);

[m,n]=size(S);

if strcmp(flag,'cubic')
    Phi= R.^3;
%     Rg = gpuArray(R);
%     Phig= Rg.^3;
%     Phi= gather(Phig);
elseif strcmp(flag,'TPS')
    R(R==0)=1;
    Phi=R.^2.*log(R);
elseif strcmp(flag,'linear')
    Phi=R;
elseif strcmp(flag,'Gaussian')
    Sd=real(sqrt(S.^2*ones(size(S'))+ones(size(S))*(S').^2-2*S*(S')));    % For Gaussian basis function
    rho=max(max(Sd))/(n*m)^(1/n);     % For Gaussian basis function
    Phi=exp(-R.^2/rho^2);
elseif strcmp(flag,'multiquad')
    Sd=real(sqrt(S.^2*ones(size(S'))+ones(size(S))*(S').^2-2*S*(S')));    % For Gaussian basis function
    rho=max(max(Sd))/(n*m)^(1/n);     % For Gaussian basis function
    Phi=sqrt((R.^2 + rho^2).^3);
elseif strcmp(flag,'invmultiquad')
    Sd=real(sqrt(S.^2*ones(size(S'))+ones(size(S))*(S').^2-2*S*(S')));    % For Gaussian basis function
    rho=max(max(Sd))/(n*m)^(1/n);     % For Gaussian basis function
    Phi=1./sqrt(R.^2 + rho^2);
end
    
Yest1=Phi*lambda; %first part of response surface
Yest2=[X,ones(mX,1)]*gamma; %optional polynomial tail
Yest=Yest1+Yest2; %predicted function value

end %functions
