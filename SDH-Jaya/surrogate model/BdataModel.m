function [BModelPoint,fd]=BdataModel(Uxdata,Ufdata,Txdata,Tfdata,xdata,fdata,mini,maxi,popsize,MaxSstag,HPop,var,gens,maxgens)
[RBFRMSE,srgtRBF] = PredictRBF(Uxdata,Ufdata,Txdata,Tfdata);
[PRRMSE,srgtPR] = PredictPR(Uxdata,Ufdata,Txdata,Tfdata);
w=[];
e=[RBFRMSE,PRRMSE];
for i=1:2
    w(i)=0.5-(e(i)/(2*(e(1)+e(2))));
end
mixfx=[fdata',xdata];
mixfx=sortrows(mixfx,1);
xdatasize=size(mixfx,1);
xpop=mixfx(1:popsize,2:end);
fpop=mixfx(1:popsize,1)'; 
amin=min(fpop);
bmin=min(fpop);
Sstag=0;
fd=0;
while(Sstag<MaxSstag)
    xnew=[];
    if rand()<(cos((gens/maxgens)*pi)+1)*0.5
        HPop=xpop;
    end
    xnew=popupdate(xpop,fpop);
    if rand()<(cos((gens/maxgens)*pi*2)+1)*0.5
        xnew=Crossover(popsize,var,xpop,HPop);
    end
    xnew=range(xnew,mini,maxi);
    YRBF = srgtsRBFEvaluate(xnew,srgtRBF);
    YPR = srgtsPRSPredictor(xnew,Uxdata,srgtPR);
    fnew=[];
    for i=1:size(xpop,1)
        fnew(i)=w(1)*YRBF(i)+w(2)*YPR(i);
    end
    for i=1:popsize
        if(fnew(i)<fpop(i))
            xpop(i,:)=xnew(i,:);
            fpop(i)=fnew(i);
        end
    end
    bmin=amin;
    amin=min(fpop);
    if bmin==amin
        Sstag=Sstag+1;
    else
        Sstag=0;
    end
end
[t,tindex]=min(fpop);
BModelPoint=xpop(tindex,:);
[t,tindex]=min(fdata);
best=xdata(tindex,:);
if min(fpop)<min(fdata)
    if best~=BModelPoint
        fd=1;
    end
end
end