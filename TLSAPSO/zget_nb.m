function [keys_out, values_out] = zget_nb(dict,key,nb_range)

% 1. 调用方法：
%   [keys_nb, values_nb] = zget_all(hash_table,key,nb_range)
% 2. 参数
%   ditc：哈希表；key：位置向量；nb_range：对应key各维的邻域范围向量
%         要求key、nb_range以及哈希表中key向量中的元素数量相等，否则调用失败，显示提示信息
% 3. 返回值
%   keys_out：cell类型，其中每个元素是一个位置向量
%             使用方法：keys_out{i}取出第i个位置向量
%                     keys_out{i}(j)取出第i个位置向量的第j个数据
%   values_out：cell类型，其中每个元素是一个适应值(double)
%               使用方法：values_out{i}取出第i个适应值

keys_nb = cell(1,1);
values_nb = cell(1,1);

len_key = size(key, 2);
len_range = size(nb_range, 2);

if len_key ~= len_range
    display('ERROR: zget_nb calling! -- Length of "key" and "nb_range" not equal!');
    keys_out = keys_nb;
values_out = values_nb;
    return;
end

[keys, values] = zget_all(dict);

len_keys = numel(keys);
if len_keys <= 0
    keys_out = keys_nb;
    values_out = values_nb;
    return;
end

len_keys_in_cell = numel(str2num(str2mat(keys(1))));
if len_keys_in_cell ~= len_range
    display('ERROR: zget_nb calling! -- Length of "key" is not equals to keys in hash!');
    keys_out = keys_nb;
values_out = values_nb;
    return;
end

nb_cur_index = 1;
for i = 1:len_keys
    key_nb = str2num(str2mat(keys(i)));
    for j = 1:len_key
        diff = abs(key_nb(j) - key(j));
        if diff > nb_range(j)
            break;
        end
    end
    if j >= len_key
        keys_nb{nb_cur_index} = key_nb;
        values_nb{nb_cur_index} = values(i);
        nb_cur_index = nb_cur_index + 1;
    end
end

keys_out = keys_nb;
values_out = values_nb;