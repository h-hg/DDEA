disttoarchive=zeros(size(archive,1),1);
for i1=1:size(archive,1)
    disttoarchive(i1,1)=0;
    for t=1:dimension
        disttoarchive(i1,1)=disttoarchive(i1,1)+(current_position(i,t)-archive(i1,t+3))^2;
    end
end
min(sqrt(disttoarchive))