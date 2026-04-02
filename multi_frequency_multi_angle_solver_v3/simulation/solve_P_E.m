function [P, E, Escat] = solve_P_E(k, afun, Chi, precon, Cpp, Einc, tol, maxit)   
% Solve polarization strength and electric fields in the simulation region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: wavenumber, 2*pi/wavelength, size (N, 1)
% afun: function of polarization currents P, P has a size of (nx*ny,1), Einc = Etot - Escat, P =
% k^2*chi*E, Escat = G*P, G is the Green's function convolution operator
% Chi: -1*scattering potential, size (ny, nx)
% precon: function that applys the inverse of (I+BG) as a preconditioner,
% ff is k^2*chi*Einc, (I+BG)^-1*ff = P
% Cpp: Green's function convolution operator without k^2 prefactor (k^2 is
% contained in P), size (2*ny-1, 2*nx-1)
% Einc: incident electric fields, size (ny, nx)
% tol: tolerance for bicgstab
% maxit: maximal iteration number for bicgstab

% Outputs:
% P: polarization currents in the simulation region,P has no k^2 prefactor,
% it has the unit of dipole amplitude, size (ny, nx)
% E: total electric field in the simulation region, size (ny, nx)
% Escat: scattered electric field in the simulation region, size (ny, nx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The output P is dipole amplitude with no k^2 prefactor
y = -reshape(Chi .* Einc,[],1);
tic; P = bicgstab(afun, y, tol, maxit,precon); toc % solve P for a given set of k, Einc and chi
escat = BTTB_matvec(Cpp, P); % scattered fields
Escat = reshape(escat, size(Einc));
E = Escat + Einc; % total fields
P = reshape(P, size(Einc))/k^2; % dipole amplitude, no k^2 prefactor
end