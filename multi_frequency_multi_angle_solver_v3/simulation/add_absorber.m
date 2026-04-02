function [chi,abs_patch_x,abs_patch_y] = add_absorber(k,target_R,absl,n_core,xSim,ySim,x1,x2,y1,y2,h,chi)
% add absorbers on the two sides of the waveguide to suppress reflection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: wavenumber, scalar
% target_R: scalar, this determines max(Im chi), smaller values for larger
% absorption
% absl: length of absorber
% n_core: waveguide core refractive index
% xSim, ySim: simulation grid mesh, size (ny, nx) 
% x1,x2,y1,y2: coordinates of waveguide corners
% h: grid size
% chi: susceptibility matrix in the simulation region,size (ny,nx)

% Outputs:
% chi: suceptibility maxtirx in the simulation region with added absorbers,
% chi has imaginary parts at the two sides of the waveguide
% abs_patch_x, abs_patch_y: coordinates to draw patches

% Notes:
% When there are multiple wavelengths, set the absorber for the largest
% wavelength and use it for all wavelengths
% The absorbers work by gradually turning on electric conductivities. Here 
% I use a quadratic function. 
% Align the positions with the simulation grid and make sure the number of 
% grid points is correct,i.e., number of points*grid size = structure size
% Ref:Steven G. Johnson, Optics Express, 2008, vol.16. no. 15.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma0 = -log(target_R)/(4*absl*n_core*1/3); % 1/3 for quadratic distribution
% 4 corners of the left absorber, rounded to align with the simulation grids
x1_left = round(x1/h)*h; 
x2_left = round((x1+absl)/h)*h; 
y1_left = round(y1/h)*h; 
y2_left = round(y2/h)*h; 
left_absorber_region = xSim>x1_left-h/2 & xSim<x2_left-h/2 & ySim>y1_left-h/2 & ySim<y2_left-h/2; % left absorber index, h/2 is needed to get the correct number of points
sigma_left = (xSim(left_absorber_region)-x2_left).^2/absl^2*sigma0/k; % quadratic distribution
chi(left_absorber_region) = (chi(left_absorber_region)+1).*(1+1i*sigma_left)-1; % eps = chi+1

% 4 corners of the right absorber, rounded to align with the simulation grids
x1_right = round((x2-absl)/h)*h; 
x2_right = round(x2/h)*h; 
y1_right = y1_left; 
y2_right = y2_left; 
right_absorber_region = xSim>x1_right-h/2 & xSim<x2_right-h/2 & ySim>y1_right-h/2 & ySim<y2_right-h/2; % right absorber index
sigma_right = (xSim(right_absorber_region)-x1_right).^2/absl^2*sigma0/k;
chi(right_absorber_region) = (chi(right_absorber_region)+1).*(1+1i*sigma_right)-1; % eps = chi+1

% coordinates of absorbers to draw patches
abs_patch_x = [x1_left,x2_right;x2_left,x1_right;x2_left,x1_right;x1_left,x2_right];
abs_patch_y = [y1_left,y1_right;y1_left,y1_right;y2_left,y2_right;y2_left,y2_right];
end