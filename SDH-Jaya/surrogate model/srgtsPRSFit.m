function srgtSRGT = srgtsPRSFit(srgtOPT)
%Function srgtsPRSFit fits the specified polynomial response surface model.
% 
%    srgtSRGT = srgtsPRSFit(srgtOPT)
% 
%srgtSRGT is the surrogate structure that contains the following fields:
%* NbPoints       : number of points in the data set.
%* NbVariables    : number of input variables.
%* PRS_Degree     : degree of the full polynomial surface that fits the data.
%* PRS_Beta       : the vector of coefficients of the polynomial response
%                   surface.
%* PRS_SE         : PRS standard error.
%* PRS_RemovedIdx : vector of indexes of those elements removed during the
%                   stepwise regression.
%
%Example:
%     % basic information about the problem
%     myFN = @cos;  % this could be any user-defined function
%     designspace = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     % create DOE
%     npoints = 5;
%     X = linspace(designspace(1), designspace(2), npoints)';
%
%     % evaluate analysis function at X points
%     Y = feval(myFN, X);
%
%     % fit surrogate models
%     options = srgtsPRSSetOptions(X, Y);
% 
%     surrogate = srgtsPRSFit(options)
% 
%     surrogate = 
% 
%           NbPoints: 5
%        NbVariables: 1
%         PRS_Degree: 2
%           PRS_Beta: [3x1 double]
%             PRS_SE: 0.3381
%     PRS_RemovedIdx: []

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate Gramian matrix X
[srgtSRGT.NbPoints srgtSRGT.NbVariables] = size(srgtOPT.P);
srgtSRGT.PRS_Degree = srgtOPT.PRS_Degree;
X = srgtsPRSCreateGramianMatrix(srgtOPT.P, srgtSRGT.NbVariables, srgtOPT.PRS_Degree);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute beta coefficients and statistics
switch srgtOPT.PRS_Regression
    case 'Full'
        [beta, BINT, R, RINT, STATS] = regress(srgtOPT.T, X);
        SE = sqrt(STATS(4));
        PRS_RemovedIdx = [];
        
    case 'StepwiseSRGTS'
        y = srgtOPT.T;
        % construct matrix X
        Xinv = pinv(X'*X);
        
        % compute beta coefficients of the full model
        beta = Xinv*X'*y;
        NbCoeff = length(beta);
       
        % get the t-statistic for the selected model
        tstatistic = (beta./sqrt(((y'*y - (beta'*X')*y)/(srgtSRGT.NbPoints - NbCoeff))*diag(Xinv)))';
        
        [mn_t, idx] = min(abs(tstatistic));
        ctr = 1;
        PRS_RemovedIdx    = [];
        while (abs(mn_t) < 1)    
            PRS_RemovedIdx(ctr) = idx;
            for c1 = 1:ctr-1
                PRS_RemovedIdx(ctr) = PRS_RemovedIdx(ctr) +(PRS_RemovedIdx(ctr) >= PRS_RemovedIdx(c1));
            end
            PRS_RemovedIdx = sort(PRS_RemovedIdx);
            ctr = ctr + 1;
            X(:,idx) = [];
            % construct matrix Xinv
            Xinv = pinv(X'*X);
            
            % compute beta coefficients of the selected model
            beta = Xinv*X'*y;
            NbCoeff = length(beta);
            
            % get the t-statistic for the selected model
            tstatistic = (beta./sqrt(((y'*y - (beta'*X')*y)/(srgtSRGT.NbPoints - NbCoeff))*diag(Xinv)))';
            
            [mn_t, idx] = min(abs(tstatistic));
        end
        SE = sqrt( (y'*y - (beta'*X')*y)/(srgtSRGT.NbPoints - NbCoeff) );
        
    case 'StepwiseMATLAB'
        [beta, se, PVAL, InModel, STATS] = stepwisefit(X, srgtOPT.T, 'display', 'off');
        PRS_RemovedIdx = find((horzcat((STATS.intercept ~= 0), InModel)') == 0);
        beta = vertcat(STATS.intercept, beta); beta(PRS_RemovedIdx) = [];
        SE   = (STATS.rmse); % standard error
        
    case 'ZeroIntercept'
        X(:,1) = [];
        [beta, STATS] = regress(srgtOPT.T, X);
        SE = sqrt(STATS(4));
        PRS_RemovedIdx = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign srgtSRGT outputs
srgtSRGT.PRS_Beta       = beta;
srgtSRGT.PRS_SE         = SE;
srgtSRGT.PRS_RemovedIdx = PRS_RemovedIdx;

return
