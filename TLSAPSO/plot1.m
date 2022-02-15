clc
clear all
global dimension
global popsize
global initial_flag
dimension=2;
popsize=60;
xmax=100;
xmin=-100;
current_position=zeros(popsize,dimension);
current_fitness=zeros(popsize,1);
global_flag=1;

dxmin=xmin;
dxmax=dxmin+(xmax-xmin)/popsize;

for i=1:popsize
    for t=1:dimension
        current_position(i,t)=rand*(dxmax-dxmin)+dxmin;
%         current_position(i,t)=rand*(xmax-xmin)+xmin;
    end
    dxmin=dxmax;
    dxmax=dxmin+(xmax-xmin)/(popsize);
    temp_position=current_position(i,:);
    initial_flag=0;
    current_fitness(i,1)=benchmark_func(temp_position,global_flag);
    currpos_arch(i,1)=i;
    currpos_arch(i,2)=0;
    currpos_arch(i,3)=current_fitness(i,1);
    for t=1:dimension
        currpos_arch(i,t+3)=current_position(i,t);
    end
end

% [x,y]=(current_position(:,1)',current_position(:,2)');
% reshape(current_fitness,size(current_fitness',1))
% z=ones(size(current_fitness))*current_fitness;
% surf(x,y,z);


plot3(current_position(:,1)',current_position(:,2)',current_fitness(:,1)','r-');
hold on;

xmax_t=zeros(dimension,1);
xmin_t=zeros(dimension,1);
for t=1:dimension
    xmax_t(t,1)=max(current_position(:,t));
    xmin_t(t,1)=min(current_position(:,t));
    nbthreshold(1,t)=abs(xmax_t(t,1)-xmin_t(t,1));
    xmax_t(t,1)=xmax_t(t,1)+0.25*(xmax_t(t,1)-xmin_t(t,1));%(xmax-xmin)/popsize;%
    xmax_t(t,1)=min(xmax,xmax_t(t,1));
    xmin_t(t,1)=xmin_t(t,1)-0.25*(xmax_t(t,1)-xmin_t(t,1));%(xmax-xmin)/popsize;%
    xmin_t(t,1)=max(xmin,xmin_t(t,1));
end


        
curr_index=1;
for i1=1:size(currpos_arch,1)
    input_curr=1;
    for t=1:dimension
        if currpos_arch(i1,t+3)<xmin_t(t,1) || currpos_arch(i1,t+3)>xmax_t(t,1)
            input_curr=0;
            break;
        end
    end
    if input_curr==1
        for t=1:dimension
            curr_pos(curr_index,t)=currpos_arch(i1,t+3);
        end
        curr_fitness(curr_index,1)=currpos_arch(i1,3);
        curr_index=curr_index+1;
    end
end
max_pos=zeros(dimension,1);
min_pos=zeros(dimension,1);
for t=1:dimension
    max_pos(t,1)=max(curr_pos(:,t));
    min_pos(t,1)=min(curr_pos(:,t));
end
maxmin_temp=abs(max_pos-min_pos);
spread1=min(maxmin_temp(find(maxmin_temp>0)));
sample_size=size(currpos_arch,1);
add_node=20;
timenet=newrb(curr_pos',curr_fitness',0.1,spread1,add_node,5);
for i=1:popsize
current_fitness1(i,1)=sim(timenet,current_position(i,:)');
end
plot3(current_position(:,1)',current_position(:,2)',current_fitness1','o');
