
function out=fitness(in,func_id) %func_id represents which function will be used
out=0;
f1=0;

if func_id == 1
    f2=1;               %f8  F5 Griewank 50dimension
    for i=1:50
        f1=f1+in(:,i)^2;
        f2=f2*cos(in(:,i)/sqrt(i));
    end
    out=1/4000*f1-f2+1;
else if func_id == 2
        f2=0;               % F4 Ackley function 50dimension
        for i=1:50
            f1=f1+in(:,i)^2;
            f2=f2+cos(2*pi*in(:,i));
        end
        out=-20*exp(-0.2*sqrt(1/50*f1))-exp(1/50*f2)+20+exp(1);
    else if func_id == 3
            for j=1:49        % F2Rosenbrock 50dimension
                f1=100*(in(:,j+1)-in(:,j)^2)^2+(in(:,j)-1)^2;
                out=out+f1;
            end
        else if func_id == 4
                for i=1:50 %Ellipsoid problem 50d
                    f1 = i * in(:,i)^2;
                    out = out + f1;
                end
            else if func_id == 5
                    f2=1;               %f8  F5 Griewank 100dimension
                    for i=1:100
                        f1=f1+in(:,i)^2;
                        f2=f2*cos(in(:,i)/sqrt(i));
                    end
                    out=1/4000*f1-f2+1;
                else if func_id == 6
                        f2=0;               % F4 Ackley function 100dimension
                        for i=1:100
                            f1=f1+in(:,i)^2;
                            f2=f2+cos(2*pi*in(:,i));
                        end
                        out=-20*exp(-0.2*sqrt(1/100*f1))-exp(1/100*f2)+20+exp(1);
                    else if func_id == 7
                            for j=1:99        % F2Rodenbrock 100dimension
                                f1=100*(in(:,j+1)-in(:,j)^2)^2+(in(:,j)-1)^2;
                                out=out+f1;
                            end
                        else if func_id == 8
                                for i=1:100 %Ellipsoid problem 100d
                                    f1 = i * in(:,i)^2;
                                    out = out + f1;
                                end
                            else if func_id == 13
                                    f2=1;               %f8  F5 Griewank 200dimension
                                    for i=1:200
                                        f1=f1+in(:,i)^2;
                                        f2=f2*cos(in(:,i)/sqrt(i));
                                    end
                                    out=1/4000*f1-f2+1;
                                else if func_id == 14
                                        f2=0;               % F4 Ackley function 200dimension
                                        for i=1:200
                                            f1=f1+in(:,i)^2;
                                            f2=f2+cos(2*pi*in(:,i));
                                        end
                                        out=-20*exp(-0.2*sqrt(1/200*f1))-exp(1/200*f2)+20+exp(1);
                                    else if func_id == 15
                                            for j=1:199        % F2Rodenbrock 200dimension
                                                f1=100*(in(:,j+1)-in(:,j)^2)^2+(in(:,j)-1)^2;
                                                out=out+f1;
                                            end
                                        else if func_id == 16
                                                for i=1:200 %Ellipsoid problem 200d
                                                    f1 = i * in(:,i)^2;
                                                    out = out + f1;
                                                end
                                            else if func_id == 17
                                                    f2=1;               %f8  F5 Griewank 500dimension
                                                    for i=1:500
                                                        f1=f1+in(:,i)^2;
                                                        f2=f2*cos(in(:,i)/sqrt(i));
                                                    end
                                                    out=1/4000*f1-f2+1;
                                                else if func_id == 18
                                                        f2=0;               % F4 Ackley function 500dimension
                                                        for i=1:500
                                                            f1=f1+in(:,i)^2;
                                                            f2=f2+cos(2*pi*in(:,i));
                                                        end
                                                        out=-20*exp(-0.2*sqrt(1/500*f1))-exp(1/500*f2)+20+exp(1);
                                                    else if func_id == 19
                                                            for j=1:499        % F2Rodenbrock 500dimension
                                                                f1=100*(in(:,j+1)-in(:,j)^2)^2+(in(:,j)-1)^2;
                                                                out=out+f1;
                                                            end
                                                        else if func_id == 20
                                                                for i=1:500 %Ellipsoid problem 500d
                                                                    f1 = i * in(:,i)^2;
                                                                    out = out + f1;
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


% for i=1:20 %Ellipsoid problem 30d
%     f1 = i * in(:,i)^2;
%     out = out + f1;
% end

% for i=1:30 %Ellipsoid problem 30d
%     f1 = i * in(:,i)^2;
%     out = out + f1;
% end

% for i=1:50 %Ellipsoid problem 50d
%     f1 = i * in(:,i)^2;
%     out = out + f1;
% end

% for i=1:100 %Ellipsoid problem 100d
%     f1 = i * in(:,i)^2;
%     out = out + f1;
% end

% for i=1:30
%     f1=in(:,i)^2;
%     out=out+f1;
% end

% for i=1:30
%     f1=abs(in(:,i))^(i+1);    %f1
%     out=out+f1;
% end

%out=100*(in(:,2)^2-in(:,1))^2+(1-in(:,1))^2;  %f2
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

%  out=4*in(:,1)^2-2.1*in(:,1)^4+1/3*in(:,1)^6+in(:,1)*in(:,2)-4*in(:,2)^2+4*in(:,2)^4+2;  %f4  F7six_hump camel_back

%   out=(1+(in(:,1)+in(:,2)+1)^2*(19-14*in(:,1)+3*in(:,1)^2-14*in(:,2)+6*in(:,1)*in(:,2)+3*in(:,2)^2))*(30+(2*in(:,1)-3*in(:,2))^2*(18-32*in(:,1)+12*in(:,1)^2+48*in(:,2)-36*in(:,1)*in(:,2)+27*in(:,2)^2));
%f5  F8 Goldstein_price

%  for i=1:30
%     f1=f1+in(:,i)*sin(sqrt(abs(in(:,i))));
%  end
%  out=12569.5-f1;   %f6
%
%  for i=1:30
%      out=out+(in(:,i)^2-10*cos(2*pi*in(:,i))+10);   %f7 F3 Rastrigin
%  end

%   f2=1;               %f8  F5 Griewank 20dimension
%   for i=1:20
%       f1=f1+in(:,i)^2;
%       f2=f2*cos(in(:,i)/sqrt(i));
%   end
%   out=1/4000*f1-f2+1;

%   f2=1;               %f8  F5 Griewank
%   for i=1:30
%       f1=f1+in(:,i)^2;
%       f2=f2*cos(in(:,i)/sqrt(i));
%   end
%   out=1/4000*f1-f2+1;


%   f2=1;               %f8  F5 Griewank 50dimension
%   for i=1:50
%       f1=f1+in(:,i)^2;
%       f2=f2*cos(in(:,i)/sqrt(i));
%   end
%   out=1/4000*f1-f2+1;

%   f2=1;               %f8  F5 Griewank 100dimension
%   for i=1:100
%       f1=f1+in(:,i)^2;
%       f2=f2*cos(in(:,i)/sqrt(i));
%   end
%   out=1/4000*f1-f2+1;
% % %
% f2=0;
% for i=1:30      % F1schwefel1.2
%     for j=1:i
%         f1=in(:,j)+f1;
%     end
%     f2=f1^2;
%     out=out+f2;
% end
% % %
% for j=1:29        % F2Rodenbrock 30dimension
%     f1=100*(in(:,j+1)-in(:,j)^2)^2+(in(:,j)-1)^2;
%     out=out+f1;
% end
% for j=1:49        % F2Rodenbrock 50dimension
%     f1=100*(in(:,j+1)-in(:,j)^2)^2+(in(:,j)-1)^2;
%     out=out+f1;
% end

% for j=1:99        % F2Rodenbrock 100dimension
%     f1=100*(in(:,j+1)-in(:,j)^2)^2+(in(:,j)-1)^2;
%     out=out+f1;
% end
%
% f2=0;               % F4 Ackley function 50dimension
% for i=1:50
%     f1=f1+in(:,i)^2;
%     f2=f2+cos(2*pi*in(:,i));
% end
% out=-20*exp(-0.2*sqrt(1/50*f1))-exp(1/50*f2)+20+exp(1);

% f2=0;               % F4 Ackley function 100dimension
% for i=1:100
%     f1=f1+in(:,i)^2;
%     f2=f2+cos(2*pi*in(:,i));
% end
% out=-20*exp(-0.2*sqrt(1/100*f1))-exp(1/100*f2)+20+exp(1);

% f2=0;               % F4 Ackley function 30dimension
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


% sum2=0; %weierstrass function
% for i=1:30
%     sum1=0;
%     for k=0:20
%         sum1=sum1+0.5^k*cos(2*3.14159*3^k*(in(:,i)+0.5));
%     end
%     sum2=sum2+sum1;
% end
% sum3=0;
% for k=0:20
%     sum3=sum3+0.5^k*cos(2*3.14159*3^k*0.5);
% end
% out=sum2-30*sum3;