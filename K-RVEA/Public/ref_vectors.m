function vectors = ref_vectors(M,Problem)
[~,~,p1,p2] = P_settings('K-RVEA',Problem,M);
[~,Vs] = F_weight(p1,p2,M);
Vs(Vs==0) = 0.000000001;
for i = 1:size(Vs,1)
    Vs(i,:) = Vs(i,:)./norm(Vs(i,:));
end;
vectors = Vs;