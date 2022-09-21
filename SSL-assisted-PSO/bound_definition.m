function [xmin xmax vmin vmax] = bound_definition(flag)
if (flag == 1) || (flag == 5)
    %Griewank
    xmin = -600;
    xmax = 600;
    vmin = -600;
    vmax = 600;
else if (flag == 2) || (flag == 6)
        %Ackley
        xmin = -32.768;
        xmax = 32.768;
        vmin = -32.768;
        vmax = 32.768;
    else if (flag == 3) || (flag == 7)
            %Rosenbrock
            xmin = -2.048;
            xmax = 2.048;
            vmin = -2.048;
            vmax = 2.048;
        else if (flag == 4) || (flag == 8) || (flag==21)
                %Ellipsoid
                xmin = -5.12;
                xmax = 5.12;
                vmin = -5.12;
                vmax = 5.12;
            else if flag == 9 || flag == 10 || flag == 11 || flag == 12 || flag == 13 || flag == 14
                    xmin = -5;
                    xmax = 5;
                    vmin = -5;
                    vmax = 5;
                end
            end
        end
    end 
end
end