function [ParentCobj,ParentCdec,ParentMobj,ParentMdec] = MatingSelection(CAobj,CAdec,DAobj,DAdec,N)
% The mating selection of Two_Arch2

    CAParent1 = randi(size(CAobj,1),1,ceil(N/2));
    CAParent2 = randi(size(CAobj,1),1,ceil(N/2));
    Dominate  = any(CAobj(CAParent1,:)<CAobj(CAParent2,:),2) - any(CAobj(CAParent1,:)>CAobj(CAParent2,:),2);  
    ParentCobj   = [CAobj([CAParent1(Dominate==1),CAParent2(Dominate~=1)],:);...
                 DAobj(randi(size(DAobj,1),1,ceil(N/2)),:)];
    ParentCdec = [CAdec([CAParent1(Dominate==1),CAParent2(Dominate~=1)],:);...
                 DAdec(randi(size(DAdec,1),1,ceil(N/2)),:)];
    ParentMobj   = CAobj(randi(size(CAobj,1),1,N),:);
    ParentMdec   = CAdec(randi(size(CAdec,1),1,N),:);
end