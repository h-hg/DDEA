function [keys_out, values_out] = zget_all(dict)

% 1. 调用方法：
%   [keys_out, values_out] = zget_all(hash_table)
% 2. 参数
%   ditc：哈希表
% 3. 返回值
%   keys_out：cell类型，其中每个元素是一个位置向量
%             使用方法：str2num(str2mat(keys_out(i)))取出第i个位置向量
%   values_out：向量类型，其中每个元素是一个适应值(double)
%               使用方法：values_out(i)取出第i个适应值

keys = cell(dict.size,1);
values = zeros(dict.size,1);

e = dict.keys();
i = 1; 
while e.hasMoreElements()
    keys{i} = e.nextElement();
    values(i) = dict.get(keys{i});
    i = i + 1;
end

keys_out = keys;
values_out = values;