function [ f ] = distance( Nap,ModelX )
%DISTANCE Summary of this function goes here
%   Detailed explanation goes here
D=[];
for i=1:size(ModelX,1)
    d=norm(Nap-ModelX(i,:));
    D=[D;d];
end
f=min(D);

end

