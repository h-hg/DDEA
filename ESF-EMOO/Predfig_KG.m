function [] = Predfig_KG(Global)

D = Global.D;
N = 11*D-1;
PopDec = lhsamp(N,D);
Population = INDIVIDUAL(repmat(Global.upper-Global.lower,N,1).*PopDec+repmat(Global.lower,N,1));
TDec = Population.decs; TObj = Population.objs;
[npoints,nvariables] = size(TDec);
Theta = 5.*ones(Global.M,Global.D);
Theta = (npoints^(-1/nvariables))*ones(Global.M, nvariables);

for i = 1 : Global.M
    dmodel     = dacefit(TDec(:,:),TObj(:,i),'regpoly0','corrgauss',Theta(i,:),1e-5.*ones(1,D),100.*ones(1,D));
    KG{i}   = dmodel;
    Theta(i,:) = dmodel.theta;
%     dmodel     = dacefit(TDec(:,:),TObj(:,i),'regpoly0','corrgauss',Theta(i,:));
%     KG{i}   = dmodel;
%     Theta(i,:) = dmodel.theta;
%     srgtOPTKRG  = srgtsKRGSetOptions(TDec, TObj(:,i));
%     KG{i} = srgtsKRGFit(srgtOPTKRG);
end


N2 = 50;
PopDec = lhsamp(N2,D);
Population2 = INDIVIDUAL(repmat(Global.upper-Global.lower,N2,1).*PopDec+repmat(Global.lower,N2,1));
Off2 = Population2.decs; obj = Population2.objs;
for i = 1: size(Off2,1)
    for j = 1 : Global.M
        [y(i,j),~,MSE(i,j)] = predictor(Off2(i,:),KG{j});
%         y(i,j) = srgtsKRGEvaluate(Off2(i,:),KG{j});
    end
end

createfigure3D(y,'Prediction','KG_SAEA');
% main('-algorithm',@Predfig_KG,'-problem',@DTLZ2,'-N',100,'-D',20,'-M',3)
end


