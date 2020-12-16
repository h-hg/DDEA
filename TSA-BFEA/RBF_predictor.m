function [TestNNOut]=RBF_predictor(W2,B2,Centers,Spreads,TestSamIn)

N=size(TestSamIn,1);
TestDistance = dist(Centers',TestSamIn');
TestSpreadsMat = repmat(Spreads,1,N);
TestHiddenUnitOut = radbas(TestDistance./TestSpreadsMat);
TestNNOut = W2*TestHiddenUnitOut+B2;
TestNNOut=TestNNOut';
end