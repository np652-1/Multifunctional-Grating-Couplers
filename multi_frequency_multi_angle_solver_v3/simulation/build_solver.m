function [afun, Chi, precon, Cpp] = build_solver(k, chi, ax, ay, nx, ny, Npv, Nph, NODES, param, quad_order)
% build sovler for one set of frequency and susceptibility distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: wavenumber, 2*pi/wavelength, size (N, 1)
% chi: susepctibilty, size (ny,nx)
% ax: simulation region horizontal size
% ay: simulation region vertical size
% nx: number of points in horizontal direction
% ny: number of points in vertical direction
% Npv: total vertical panel number, see identical_grid_info.m
% Nph: total horizontal panel number, see identical_grid_info.m
% NODES: simulation nodes without scatterer information from build_nodes,
% size {N,1}
% param: solver parameters and quadrature correction coefficients, size
% {N,1}, type cell, each element is a struct, output from build_nodes
% quad_order: quadrature correction for the singularities in Green's
% function, used in build_greens_function

% Outputs:
% afun: function of polarization currents P, P has a size of (nx*ny,1), Einc = Etot - Escat, P =
% k^2*chi*E, Escat = G*P, G is the Green's function convolution operator
% Chi: -1*scattering potential, size (ny, nx)
% precon: function that applys the inverse of (I+BG) as a preconditioner,
% ff is k^2*chi*Einc, (I+BG)^-1*ff = P
% Cpp: Green's function convolution operator without k^2 prefactor (k^2 is
% contained in P), size (2*ny-1, 2*nx-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = ax / (nx - 1); % grid size
Chi = -k^2*chi; % scattering potential

% loop over leaf nodes, store Chi on each panel
NODES = get_geometry(NODES,Chi(:));

% build inverse
tic
NODES = inverse_build(NODES,k,h,param); % construct the inverse of the integral operator
tbuild_solver = toc

Cpp = build_greens_function(k, ax, ay, nx, ny, quad_order, param.D0, param.D1); % construct Green's function convolution operator
afun = @(P) P + reshape(Chi,[],1) .* BTTB_matvec(Cpp, P);
precon = @(ff) apply_inverse_solve(NODES,ff,Npv,Nph); % apply the inverse as a preconditioner
end