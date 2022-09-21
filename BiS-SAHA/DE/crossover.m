function U=crossover(X,V,CR)
global Dim
[NP,~]=size(X);

U = X;
num_rand = rand(NP, Dim);
jRand=floor(rand(NP, 1) * Dim);
for i=1:NP
    for j=1:Dim
        if (num_rand(i, j)) < CR || (j == jRand(i, 1))
            U(i,j)=V(i,j);
        end     
    end    
end

