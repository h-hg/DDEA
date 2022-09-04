%% 多项式回归模型模型
function [PRRMSE,srgtPR]=PredictPR(Uxdata,Ufdata,Txdata,Tfdata)
srgtOPT  = srgtsPRSSetOptions(Uxdata,Ufdata');
srgtPR = srgtsPRSFit(srgtOPT);
[Yhat PredVar] = srgtsPRSPredictor(Txdata,Uxdata,srgtPR);
CRITERIA = srgtsErrorAnalysis(srgtOPT,srgtPR,Tfdata',Yhat);
PRRMSE=CRITERIA.RMSE;
end 