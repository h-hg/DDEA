function  A1 = UpdataArchive(A1,New,V,mu,NI)
% Update archive

    %% Delete duplicated solutions
    All       = [A1.dec;New.dec];
    [~,index] = unique(All,'rows');
    ALL.dec       = [A1.dec;New.dec];
    ALL.obj       = [A1.obj;New.obj];
%     Total     = ALL(index);
    Total.dec = ALL.dec(index,:);
    Total.obj = ALL.obj(index,:);
    %% Select NI solutions for updating the models 
	if length(Total)>NI
        [~,active] = NoActive(New.obj,V);
        Vi         = V(setdiff(1:size(V,1),active),:);
        % Select the undeplicated solutions without re-evaluated
        % solutions
        index = ismember(Total.dec,New.dec,'rows');
        Total = Total(~index);
        % Since the number of inactive reference vectors is smaller than
        % NI-mu, we cluster the solutions instead of reference vectors
        PopObj = Total.objs;
        PopObj = PopObj - repmat(min(PopObj,[],1),length(Total),1);
        Angle  = acos(1-pdist2(PopObj,Vi,'cosine'));
        [~,associate] = min(Angle,[],2);
        Via    = Vi(unique(associate)',:);
        Next   = zeros(1,NI-mu);
        if size(Via,1) > NI-mu
            [IDX,~] = kmeans(Via,NI-mu);
            for i = unique(IDX)'
                current = find(IDX==i);
                if length(current)>1
                    best = randi(length(current),1);
                else
                    best = 1;
                end
                Next(i)  = current(best);
            end
        else
            % Cluster solutions based on objective vectors when the number
            % of active reference vectors is smaller than NI-mu
            [IDX,~] = kmeans(Total.objs,NI-mu);
                for i   = unique(IDX)'
                    current = find(IDX==i);
                    if length(current)>1
                        best = randi(length(current),1);
                    else
                        best = 1;
                    end
                    Next(i)  = current(best);
                end
            A1 = [Total(Next),New];
        end
   else 
       A1 = Total;
   end
end       