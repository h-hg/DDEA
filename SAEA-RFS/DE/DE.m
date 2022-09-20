function U=DE(X,bestX)
%变异向量用函数mutation(X,bestX,F,mutationStrategy)
%交叉向量用函数crossover(X,V,CR,crossStrategy)
%mutation
%mutationStrategy=1：DE/rand/1,
%mutationStrategy=2：DE/best/1,
%mutationStrategy=3：DE/rand-to-best/1,
%mutationStrategy=4：DE/best/2,
%mutationStrategy=5：DE/rand/2.
%crossover
%crossStrategy=1:binomial crossover
%crossStrategy=2:Exponential crossover

F=0.8;%scaling factor 缩放因子
CR=1;%crossover rate 交叉概率
mutationStrategy=2;%变异策略
crossStrategy=1;%交叉策略

% while Generation<maxIteration
     V=mutation(X,bestX,F,mutationStrategy);
     U=crossover(X,V,CR,crossStrategy);
% end
end