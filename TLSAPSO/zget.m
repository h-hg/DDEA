function ret = zget(dict,key,precise)

if(true ~= isa(dict,'java.util.Hashtable'))
    ret = false;
end

keyString = '[';
for i = 1:size(key,2)
    if i > 1
        keyString = strcat(keyString,',');
    end
    
    if precise == 0
        tmp = sprintf('%0.0f',key(i));
    elseif precise == 1
        tmp = sprintf('%0.1f',key(i));
    elseif precise == 2
        tmp = sprintf('%0.2f',key(i));
    elseif precise == 3
        tmp = sprintf('%0.3f',key(i));
    elseif precise == 4
        tmp = sprintf('%0.4f',key(i));
    elseif precise == 5
        tmp = sprintf('%0.5f',key(i));
    elseif precise == 6
        tmp = sprintf('%0.6f',key(i));
    elseif precise == 7
        tmp = sprintf('%0.7f',key(i));
    elseif precise == 8
        tmp = sprintf('%0.8f',key(i));
    elseif precise == 9
        tmp = sprintf('%0.9f',key(i));
    elseif precise == 10
        tmp = sprintf('%0.10f',key(i));
    elseif precise == 11
        tmp = sprintf('%0.11f',key(i));
    elseif precise == 12
        tmp = sprintf('%0.12f',key(i));
    elseif precise == 13
        tmp = sprintf('%0.13f',key(i));
    elseif precise == 14
        tmp = sprintf('%0.14f',key(i));
    elseif precise == 15
        tmp = sprintf('%0.15f',key(i));
    elseif precise == 16
        tmp = sprintf('%0.16f',key(i));
    else
        tmp = sprintf('%0.6f',key(i));
    end
    
    keyString = strcat(keyString,tmp);
end

keyString = strcat(keyString,']');

ret = dict.get(keyString);





