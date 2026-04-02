function [x1,x2,y1,y2,region_index,nx,ny,patch_x,patch_y] = add_region(h,region_center,region_size,xSim,ySim)
% return the corner position, logical indexing, and patch position (for
% plotting) for a region with given center and size
% align the positions with the grid and make sure the number of grid points
% is correct. number of points*grid size = structure size
x1 = round((region_center(1)-region_size(1)/2)/h)*h; 
x2 = round((region_center(1)+region_size(1)/2)/h)*h; 
y1 = round((region_center(2)-region_size(2)/2)/h)*h; 
y2 = round((region_center(2)+region_size(2)/2)/h)*h; 
region_index = xSim>x1-h/2 & xSim<x2-h/2 & ySim>y1-h/2 & ySim<y2-h/2; % index for the region  
nx = size(unique(xSim(region_index)), 1); % number of horizontal grid points in this region
ny = size(unique(ySim(region_index)), 1); % % number of vertical grid points in this region
patch_x = [x1;x2;x2;x1]; % patch position for plotting
patch_y = [y1;y1;y2;y2]; % patch position for plotting
end