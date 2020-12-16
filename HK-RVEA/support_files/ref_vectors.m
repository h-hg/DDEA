function vectors = ref_vectors(M)
p1 = [ 99 13  7  5  4  3  3  2  3];
p2 = [ 0  0  0  0  1  2  2  2  2];
p1 = p1(M-1);
p2 = p2(M-1);
if M>10
    error('Define the number of ref vectors parameters for objectives > 10');
end

[~,Vs] = F_weight(p1,p2,M);
% Vs(Vs==0) = 0.000000001;
for i = 1:size(Vs,1)
    Vs(i,:) = Vs(i,:)./norm(Vs(i,:));
end
vectors = Vs;