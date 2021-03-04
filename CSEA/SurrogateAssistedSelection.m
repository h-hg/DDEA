function Next = SurrogateAssistedSelection(Global,net,p0,p1,Ref,Input,wmax,tr)
%Surrogate-assisted selection for selecting promising solutions

% Copyright 2015-2016 Cheng He
      % Offspring population
        Next    = Reproduction(Global,Ref,Input);
        Label   = net.predict(Next);
        a = tr;b=1-tr;
        i       = 0;
        if p0<0.4||(p1<a&&p0<b)
            while i<wmax
                [~,index] = sort(Label,'descend');
                Input   = Next(index(1:length(Ref)),:);
                Next    = Reproduction(Global,Ref,Input); 
                Label   = net.predict(Next);
                i = i+size(Next,1);
            end
            Next = Next(Label>0.9,:);
        elseif p0>b&&p1<a
            Next =[];
        elseif p1>b
            while i<wmax
                [~,index] = sort(Label);
                Input   = Next(index(1:length(Ref)),:);
                Next    = Reproduction(Global,Ref,Input); 
                Label   = net.predict(Next);
                i = i+size(Next,1);
            end
            Next = Next(Label<0.1,:);
        else
            Next = [];
        end
end