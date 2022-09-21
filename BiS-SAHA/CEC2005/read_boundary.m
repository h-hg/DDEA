function [a, b] = read_boundary
global initial_flag
global func_num
if initial_flag==0
    if func_num==1 a = -100; b = 100; %[-100,100]
    elseif func_num==2 a = -100; b = 100; %[-100,100]
    elseif func_num==3 a = -100; b = 100; %[-100,100]
    elseif func_num==4 a = -100; b = 100; %[-100,100]
    elseif func_num==5 a = -100; b = 100; %[no bound],initial[-100,100];
    elseif func_num==6 a = -100; b = 100; %[-100,100]
    elseif func_num==7 a = -600; b = 600; %[-600,600]
    elseif func_num==8 a = -32; b = 32;%[-32,32]
    elseif func_num==9 a = -5; b = 5; %[-5,5]
    elseif func_num==10 a = -5; b = 5; %[-5,5]
    elseif func_num==11 a = -0.5; b = 0.5; %[-0.5,0.5]
    elseif func_num==12 a = -pi; b = pi; %[-pi,pi]
    elseif func_num==13 a = -3; b = 1; %[-3,1] 
    elseif func_num==14 a = -100; b = 100; %[-100,100]
    elseif func_num==15 a = -5; b = 5;%[-5,5]
    elseif func_num==16 a = -5; b = 5; %[-5,5]
    elseif func_num==17 a = -5; b = 5; %[-5,5]        
    elseif func_num==18 a = -5; b = 5;%[-5,5]    
    elseif func_num==19 a = -5; b = 5;%[-5,5]    
    elseif func_num==20 a = -5; b = 5; %[-5,5]  
    elseif func_num==21 a = -5; b = 5; %[-5,5]    
    elseif func_num==22 a = -5; b = 5; %[-5,5]  
    elseif func_num==23 a = -5; b = 5; %[-5,5]   
    elseif func_num==24 a = -5; b = 5; %[-5,5]  
    elseif func_num==25 a = -5; b = 5; %[-5,5]  
    elseif func_num==101 a = -5.12; b = 5.12; %[-5.12,5.12]  
    elseif func_num==102 a = -2.048; b = 2.048; %[-2.048,2.048]  
    elseif func_num==103 a = -32.768; b = 32.768; %[-32.768,32.768]  
    elseif func_num==104 a = -600; b = 600; %[-600,600] 
    elseif func_num==105 a = -5; b = 5; %[-5,5]  
    end
end