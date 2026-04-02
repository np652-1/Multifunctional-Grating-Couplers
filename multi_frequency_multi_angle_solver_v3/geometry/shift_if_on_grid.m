function [shifted_x]=shift_if_on_grid(x,h)
% This function shift the point by half a grid size if the point is on grid
% h: grid size
% This is done to avoid singularities in the Green's function
if mod(x,h)<1e-6 % if on grid, shift by h/2 to avoid singularities
   shifted_x = x-h/2;
else
   shifted_x = x;
end
end
