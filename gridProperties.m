function [dx, x_num] = gridProperties(x)

x_num = numel(x);   
dx = (x(end) - x(1))/(x_num - 1);

end