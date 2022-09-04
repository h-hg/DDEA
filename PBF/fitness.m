function Population = fitness (PopDec,problem)

Population.dec = PopDec;

switch problem.problem
    case 'DTLZ1'
        g              = 100*(problem.D-problem.M+1+sum((PopDec(:,problem.M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,problem.M:end)-0.5)),2));
        Population.obj = 0.5*repmat(1+g,1,problem.M).*fliplr(cumprod([ones(size(PopDec,1),1),PopDec(:,1:problem.M-1)],2)).*[ones(size(PopDec,1),1),1-PopDec(:,problem.M-1:-1:1)];
        
end


end