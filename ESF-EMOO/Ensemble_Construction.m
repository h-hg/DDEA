function [Ensemble,r,THETA] = Ensemble_Construction(Dec,Obj,THETA,dr,Global)

[N,D]= size(Dec);
r1 = sort(randperm(D,dr));  %20D-12
r2 = sort(randperm(D,dr));
r3 = sort(randperm(D,D));
% r1 = r3;  r2 = r3;
r = {r1;r2;r3};
Theta1 = THETA.Theta1;  Theta2 = THETA.Theta2;  Theta3 = THETA.Theta3;
% Theta1 = 5.*ones(Global.M,Global.D);  Theta2 = 5.*ones(Global.M,Global.D);
Model1 = cell(1,3);  Model2 = cell(1,3);  Model3 = cell(1,3);

[f,~] = NDSort(Obj,1);
fin = find(f==1);
rin = find(f~=1);
FObj = Obj(fin,:);
FDec = Dec(fin,:);
RDec = Dec(rin,:);
RObj = Obj(rin,:);
if N > 11*D-1+25
    if length(fin) < 11*D-1+25
        [ind,~]= kmeans(RObj,11*D-1+25-length(fin));
        Next   = zeros(1,11*D-1+25-length(fin));
        for i = unique(ind)'
            current = find(ind==i);
            if length(current)>1
                best = randi(length(current),1);
            else
                best = 1;
            end
            Next(i)  = current(best);
        end
        TDec = [FDec;RDec(Next,:)];
        TObj = [FObj;RObj(Next,:)];
    else
        [ind,~]= kmeans(FObj,11*D-1+25);
        Next   = zeros(1,11*D-1+25);
        for i = unique(ind)'
            current = find(ind==i);
            if length(current)>1
                best = randi(length(current),1);
            else
                best = 1;
            end
            Next(i)  = current(best);
        end
        TDec = FDec(Next,:);
        TObj = FObj(Next,:);
    end
else
    TDec = Dec;
    TObj = Obj;
end

for i = 1 : Global.M
    dmodel     = dacefit(TDec(:,r1),TObj(:,i),'regpoly0','corrgauss',Theta1(i,r1),1e-5.*ones(1,dr),100.*ones(1,dr));
    Model1{i}   = dmodel;
    Theta1(i,r1) = dmodel.theta;
    
    dmodel     = dacefit(TDec(:,r2),TObj(:,i),'regpoly0','corrgauss',Theta2(i,r2),1e-5.*ones(1,dr),100.*ones(1,dr));
    Model2{i}   = dmodel;
    Theta2(i,r2) = dmodel.theta;
    
    dmodel     = dacefit(TDec(:,r3),TObj(:,i),'regpoly0','corrgauss',Theta3(i,r3),1e-5.*ones(1,D),100.*ones(1,D));
    Model3{i}   = dmodel;
    Theta3(i,r3) = dmodel.theta;
end
Ensemble = {Model1,Model2,Model3};
THETA.Theta1 = Theta1;  THETA.Theta2 = Theta2;  THETA.Theta3 = Theta3;
end

