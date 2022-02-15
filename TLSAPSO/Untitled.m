d=[];
for i=1:popsize-1
    for j=i+1:popsize
d(i,j)=0;
for t=1:dimension
d(i,j)=d(i,j)+(current_position(i,t)-current_position(j,t))^2;
end
d(i,j)=sqrt(d(i,j));
    end
    
end

d=[];
for i=1:popsize-1
    for j=i+1:popsize
d(i,j)=0;
for t=1:dimension
d(i,j)=d(i,j)+(pbest(i,t)-pbest(j,t))^2;
end
d(i,j)=sqrt(d(i,j));
    end
    
end