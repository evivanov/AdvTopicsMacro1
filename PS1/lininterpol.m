function [f] = lininterpol(grid, values_on_grid, x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[~,cl] = min(abs(grid - x));
if cl == length(grid)
    cl = cl-1;
end
if (grid(cl)-x)*(grid(cl+1)-x) > 0
    cl = cl - 1;
end
f = values_on_grid(cl) + (values_on_grid(cl+1) - values_on_grid(cl))/(grid(cl+1)-grid(cl)) * (x - grid(cl));
end