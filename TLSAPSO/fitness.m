
function out=fitness(in)
out=0;
f1=0;

% for i=1:30
%     f1=abs(in(:,i))^(i+1);    %f1
%     out=out+f1;
% end

% out=100*(in(:,2)^2-in(:,1))^2+(1-in(:,1))^2;  %f2
% % 
% a=[-32 -16 0 16 32 -33 -32 -16 0 16 32 -33 -32 -16 0 16 32 -33 -32 -16 0 16 32 -33 -32 -16 0 16 32;-32 -32 -32 -32 -32 -16 -32 -32 -32 -32 -32 -16 -32 -32 -32 -32 -32 -16 -32 -32 -32 -32 -32 -16 -32 -32 -32 -32 -32];
% f2=0;
% for j=1:25                   %f3 
%     for i=1:2
%        f1=f1+(in(:,i)-a(i,j))^6;
%     end
%        f2=f2+1/(j+f1);
% end
% out=(1/500+f2)^-1;

 %out=4*in(:,1)^2-2.1*in(:,1)^4+1/3*in(:,1)^6+in(:,1)*in(:,2)-4*in(:,2)^2+4*in(:,2)^4+2;  %f4  F7six_hump camel_back

out=(1+(in(:,1)+in(:,2)+1)^2*(19-14*in(:,1)+3*in(:,1)^2-14*in(:,2)+6*in(:,1)*in(:,2)+3*in(:,2)^2))*(30+(2*in(:,1)-3*in(:,2))^2*(18-32*in(:,1)+12*in(:,1)^2+48*in(:,2)-36*in(:,1)*in(:,2)+27*in(:,2)^2));
% f5  F8 Goldstein_price

%  for i=1:30       
%     f1=f1+in(:,i)*sin(sqrt(abs(in(:,i))));
%  end
%  out=12569.5-f1;   %f6
%  
%  for i=1:30
%      out=out+(in(:,i)^2-10*cos(2*pi*in(:,i))+10);   %f7 F3 Rastrigin
%  end

%   f2=1;               %f8  F5 Griewank
%   for i=1:30
%       f1=f1+in(:,i)^2;
%       f2=f2*cos(in(:,i)/sqrt(i));
%   end
%   out=1/4000*f1-f2+1;

% f2=0;  
% for i=1:30      % F1schwefel1.2
%     for j=1:i
%         f1=in(:,j)+f1;
%     end
%     f2=f1^2;
%     out=out+f2;
% end
% 
% for j=1:29        % F2Rodenbrock 30dimension
%     f1=100*(in(:,j+1)-in(:,j)^2)^2+(in(:,j)-1)^2;
%     out=out+f1;
% end

% f2=0;               % F4 Ackley function
% for i=1:30
%     f1=f1+in(:,i)^2;
%     f2=f2+cos(2*pi*in(:,i));
% end
% out=-20*exp(-0.2*sqrt(1/30*f1))-exp(1/30*f2)+20+exp(1);

%  for i=1:30       
%     f1=f1+in(:,i)*sin(sqrt(abs(in(:,i))));
%  end
%  out=-f1;   %F6schwefel2.26

%  out=0.5+((sin(sqrt(in(:,1)^2+in(:,2)^2)))^2-0.5)/(1+0.001*(in(:,1)^2+in(:,2)^2))^2;  %F9schafferF6

% a=[0.1957 0.1947 0.1735 0.1600 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246;1/0.25 1/0.5 1 1/2 1/4 1/6 1/8 1/10 1/12 1/14 1/16];
% for i=1:11
%     f1=(a(1,i)-in(:,1)*(a(2,i)^2+a(2,i)*in(:,2))/(a(2,i)^2+a(2,i)*in(:,3)+in(:,4)))^2; 
%     out=out+f1;        %F10kowalik's function
% end
%  

