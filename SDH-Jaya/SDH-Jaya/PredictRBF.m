    %% 径向基网络模型
function [RBFRMSE,srgtRBF]=PredictRBF(Uxdata,Ufdata,Txdata,Tfdata)
srgtOPT=srgtsRBFSetOptions(Uxdata, Ufdata');
srgtRBF=srgtsRBFFit(srgtOPT);
Yhat=srgtsRBFEvaluate(Txdata, srgtRBF);
CRITERIA = srgtsErrorAnalysis(srgtOPT, srgtRBF, Tfdata', Yhat);
RBFRMSE=CRITERIA.RMSE;
end
