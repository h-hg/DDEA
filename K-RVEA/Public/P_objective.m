function [Output,Boundary,Coding] = P_objective(Operation,Problem,M,Input)


    k = find(~isstrprop(Problem,'digit'),1,'last');
    switch Problem(1:k)
        case 'DTLZ'
            [Output,Boundary,Coding] = P_DTLZ(Operation,Problem,M,Input);
        otherwise
            error(['Dude provide the problem']);
    end
end

function [Output,Boundary,Coding] = P_DTLZ(Operation,Problem,M,Input)
    persistent K;
    Boundary = NaN; Coding = NaN;
    switch Operation
        case 'init'
            k = find(~isstrprop(Problem,'digit'),1,'last');               
            K = (11-M+1)*ones(7,1);                
            K = K(str2double(Problem(k+1:end)));
            D = M+K-1;
            MaxValue   = ones(1,D);
            MinValue   = zeros(1,D);
            Population = rand(Input,D);
            Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
            
            Output   = Population;
            Boundary = [MaxValue;MinValue];
            Coding   = 'Real';
        
        case 'value'
            Population    = Input;
            K = (11-M+1);
            FunctionValue = zeros(size(Population,1),M);
            switch Problem
                case 'DTLZ2'
                    g = sum((Population(:,M:end)-0.5).^2,2);
                    for i = 1 : M
                        FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                        if i > 1
                            FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                        end
                    end
            end
            Output = FunctionValue;
       
        case 'true'
            switch Problem
                case {'DTLZ2'}
                    Population = T_uniform(Input,M);
                    for i = 1 : size(Population,1)
                    	Population(i,:) = Population(i,:)./norm(Population(i,:));
                    end           
            end
            Output = Population;
    end
end


function W = T_uniform(k,M)
    H = floor((k*prod(1:M-1))^(1/(M-1)));
    while nchoosek(H+M-1,M-1) >= k && H > 0
        H = H-1;
    end
    if nchoosek(H+M,M-1) <= 2*k || H == 0
        H = H+1;
    end
    k = nchoosek(H+M-1,M-1);
    Temp = nchoosek(1:H+M-1,M-1)-repmat(0:M-2,nchoosek(H+M-1,M-1),1)-1;
    W = zeros(k,M);
    W(:,1) = Temp(:,1)-0;
    for i = 2 : M-1
        W(:,i) = Temp(:,i)-Temp(:,i-1);
    end
    W(:,end) = H-Temp(:,end);
    W = W/H;
end

