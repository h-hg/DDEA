% <problem> <DTLZ>
% Benchmark MOP proposed by Deb, Thiele, Laumanns, and Zitzler

%------------------------------- Reference --------------------------------
% K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, Scalable test problems
% for evolutionary multiobjective optimization, Evolutionary multiobjective
% Optimization. Theoretical Advances and Applications, 2005, 105-145.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
        %% Calculate objective values
        function PopObj = DTLZ7M2(PopDec,m)
            M               = m;
            PopObj          = zeros(size(PopDec,1),M);
            g               = 1+9*mean(PopDec(:,M:end),2);
            PopObj(:,1:M-1) = PopDec(:,1:M-1);
            PopObj(:,M)     = (1+g).*(M-sum(PopObj(:,1:M-1)./(1+repmat(g,1,M-1)).*(1+sin(3*pi.*PopObj(:,1:M-1))),2));
        end