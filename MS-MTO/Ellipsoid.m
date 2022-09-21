function obj = Ellipsoid(var)
%Ellipsoid function
%   - var: design variable vector

    obj=0;
    for i=1:length(var)
    obj=obj+var(i)*var(i)*i;
    end
end