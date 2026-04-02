function [Cpp] = build_greens_function(k, ax, ay, nx, ny, quad_order, D0, D1)
% build Green's function convolutor for a given k and simulation region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: wavenumber, 2*pi/wavelength, size (N, 1)
% ax: simulation region horizontal size
% ay: simulation region vertical size
% nx: number of points in horizontal direction
% ny: number of points in vertical direction
% quad_order: quadrature correction for the singularities in Green's function 
% Notice that in the build_nodes code, we already have a
% param.quad_order. That one is used in inverse_build. Here we construct
% Green's function convolutor, and it can have a different quadrature
% correction order than the inverse.
% D0: 4th order correction coefficients
% D1: 6th order correction coefficients

% Outputs:
% Cpp: Green's function convolution operator without k^2 prefactor (k^2 is
% contained in P), size (2*ny-1, 2*nx-1). With Cpp, we can use FFT to run
% fast convolution calculation with polarization currents to get scattered
% fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hx = ax / (nx - 1); 
hy = ay / (ny - 1); % This looks unnecessary but the identical_grid_info is not perfect. hx and hy could have very small difference and leads to error in x and y, then Cpp. We want the length of x to be 2*nx-1, and length of y to be 2*ny-1
% create the GF convolutor
x = -ax:hx:ax; 
y = -ay:hy:ay;
[X, Y] = meshgrid(x,y);
r = sqrt(X.^2 + Y.^2);
g = 1j*k^2 / 4 * besselh(0,k*r);
g(ny,nx) = 0; % set singularity (at the origin) to zero

% high-order quadrature correction
if abs(hx-hy)>1e-16
    warning('Current high-order quadrature correction requires: hx = hy')
end
Tao = zeros(size(g));
switch quad_order
    case 2 % no correction, do nothing
    case 4 % forth order correction
        Tao(ny,nx) = D0;
    case 6 % six order correction
        Tao(ny,nx) = D0 - 2*D1;
        Tao(ny-1,nx) = 1/2*D1;
        Tao(ny+1,nx) = 1/2*D1;
        Tao(ny,nx-1) = 1/2*D1;
        Tao(ny,nx+1) = 1/2*D1;
end
g = g + 1j*k^2/4*Tao;
g = g*hx^2/k^2; % h^2 for integration, /k^2: because P contains the k^2 prefactor  
gg = ifftshift(g);
Cpp = fft2(gg); % Cpp does not have k^2 prefactor
end