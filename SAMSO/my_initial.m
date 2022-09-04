function[Xtrain Ytrain]=my_initial(Dat,Doe)
%根据样本点个数、实验设计方法、问题信息设计初始样本点
switch func2str(Doe)
    
    case 'DOEOLHS'       %最优拉丁超立方设计
        npoints=Dat.npoints;
        ndv=Dat.ndv;
        method= 'ESEA';
        LHS = DOEOLHS(npoints, ndv, method);
        LHS = ScaleVariable(LHS, [0 1]'*ones(1,Dat.ndv), Dat.designspace);
        Xtrain=LHS;
        Ytrain=feval(Dat.myFN, Xtrain);
        
    case 'DOETPLHS'       %平移拉丁超立方设计
        npoints=Dat.npoints;
        ndv=Dat.ndv;
        LHS =DOETPLHS(nPoints, ndv);
        LHS = ScaleVariable(LHS , [0 1]'*ones(1,Dat.ndv), Dat.designspace);
        Xtrain=LHS;
        Ytrain=feval(Dat.myFN, Xtrain);
    case 'DOELHS'       %拉丁超立方设计
        npoints=Dat.npoints;
        ndv=Dat.ndv;
        LHS =DOELHS(npoints, ndv,5);
        LHS = ScaleVariable(LHS , [0 1]'*ones(1,Dat.ndv), Dat.designspace);
        Xtrain=LHS;
        Ytrain=feval(Dat.myFN, Xtrain);
        
end


        
