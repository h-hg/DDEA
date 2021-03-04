function CSEA(Global)
% <algorithm> <2016>
%A Classification based surrogate-assisted MaOEA (CSEA)

% Copyright 2015-2016 Cheng He
    [k,gmax,H,Iteration] = Global.ParameterSet(6,3000,ceil(Global.D*2),800);
    %% Initalize the population by Latin hypercube sampling
    N           = 11*Global.D-1;
    PopDec      = lhsamp(N,Global.D);
    Population  = INDIVIDUAL(repmat(Global.upper-Global.lower,N,1).*PopDec+repmat(Global.lower,N,1));
    Arc         = Population;
	while Global.NotTermination(Arc)
	%% Select reference solutions and preprocess the data
        Ref     = RefSelect(Population,k);
        Input   = Population.decs;  
        Output  = GetOutput(Population.objs,Ref.objs); 
        rr      = sum(Output)/length(Output);
        tr      = min(rr,1-rr)*0.5;
        [TrainIn,TrainOut,TestIn,TestOut] = DataProcess(Input,Output);
	%% Construct and update the FNN
        net = NN(H,Iteration);
        net.train(TrainIn,TrainOut);
    %% Error rates calculation
        TestPre     = net.predict(TestIn);
        IndexGood   = TestOut==1;
        p0  = sum(abs((TestOut(IndexGood)-TestPre(IndexGood))))/sum(IndexGood);
        p1  = sum(abs((TestOut(~IndexGood)-TestPre(~IndexGood))))/sum(~IndexGood);
    %% Surrogate-assisted selection and update the population
         Next = SurrogateAssistedSelection(Global,net,p0,p1,Ref,Population.decs,gmax,tr);
        if ~isempty(Next)
            Arc = [Arc,INDIVIDUAL(Next)];
        end
        Population = RefSelect(Arc,Global.N);
	end
end

