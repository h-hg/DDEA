function [ d ] = fdistance( a,A )
%FDISTANCE 此处显示有关此函数的摘要
%   此处显示详细说明
if isempty(A)
    d=1;
else
    d=min(pdist2(a,A));
end

end

