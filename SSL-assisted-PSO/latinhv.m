%拉丁超立方体采样
function [current_position] = latinhv(popsize,dimension,xmin,xmax)

latin_position = zeros(dimension, popsize);
sample_selected = zeros(dimension, popsize);
for i = 1:dimension
    for j = 1:popsize
        latin_position(i,j) = xmin + ((xmax-xmin)/popsize) * j;
        sample_selected(i,j) = 0;
    end
end

pop_temp = popsize;

for i = 1:popsize
    if(pop_temp == 1)
        for j = 1:dimension
            initial_position(i,j) = latin_position(j,1);
        end
    else
        latin_position_temp = zeros(dimension, pop_temp - 1);
        for j = 1:dimension
            selected_value = fix(rand() * pop_temp);
            while (selected_value == 0)
                selected_value = fix(rand() * pop_temp);
            end
            initial_position(i,j) = latin_position(j,selected_value);
            %            sample_selected(j,selected_value) = 1;
            a = latin_position(j,:);
            b = latin_position(j,selected_value);
            c = setdiff(a,b);
            latin_position_temp(j,:) = c;
        end
        latin_position = zeros(dimension, pop_temp);
        latin_position = latin_position_temp;
        pop_temp = pop_temp - 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%分配种群
for i=1:popsize
    for t=1:dimension
        current_position(i,t)=initial_position(i,t);
    end
end

end