function [ I] = StochasticRanking( R1,R2)
I=[1:size(R1,1)]';
T=[I,R1,R2];
for i=1:size(R1,1)
    for j=1:size(R1,1)-1
       if rand>0.5
           if T(j,2)>T(j+1,2)
               t=T(j,:);
               T(j,:)=T(j+1,:);
               T(j+1,:)=t;
           end
       else
           if T(j,3)>T(j+1,3)
               t=T(j,:);
               T(j,:)=T(j+1,:);
               T(j+1,:)=t;
           end
       end
    end
end
I=T(:,1);

end

