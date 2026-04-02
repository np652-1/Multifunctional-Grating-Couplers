function [Px, Einc_adj] = poynting_flux_Px_fom(k, xySim, P, hs, Ez_inc, Hy_inc, recPts, hr, coeff, norm)
% This function calculates the integral of poynting flux Px over a line monitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: wavenumber, 2*pi/wavelength, must be a scalar
% xySim: x y points in the simulation region, size (N1,2)
% P: polarization fields in the simulation region, size (ny, nx), no k^2
% prefactor, N1 = ny*nx
% hs: grid size in the source region
% Ez_inc: incident field Ez on recPts, size (N2,1)
% Hy_inc: incident field Hy on recPts, size (N2,1)
% recPts: monitor points, size (N2,2), recPts(:,1) x, recPts(:,2) y
% hr: grid size on the monitor
% coeff: +1 for positive x direction, -1 for negative x direction
% norm: normalzation factor

% Outputs:
% Px: poynting flux in x direction, scalar
% Einc_adj: adjoint incident fields, size (ny,nx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ez_scat = electric_dipole_z_field_Ez(k, recPts, xySim, hs^2, P); % Ez fields radiated by P
Hy_scat = electric_dipole_z_field_Hy(k, recPts, xySim, hs^2, P); % Hy fields radiated by P
Ez = Ez_scat+Ez_inc; % total Ez field
Hy = Hy_scat+Hy_inc; % total Hy field
Px = -coeff*1/2*sum(real(Ez.*conj(Hy)*hr))/norm; % Poynting flux in x direction

pz_adj = -coeff*1/4*conj(Hy)/norm; % electric dipole source amplitude
my_adj = coeff*1/4*conj(Ez)/norm; % the mangetic dipole has different sign than the electric dipole
Ez_p = electric_dipole_z_field_Ez(k, xySim, recPts, hr, pz_adj); % Ez fields radiated by pz_adj
Ez_m = magnetic_dipole_y_field_Ez(k, xySim, recPts, hr, my_adj); % Ez fields radiated by my_adj
Einc_adj = Ez_p + Ez_m;
Einc_adj = reshape(Einc_adj, size(P));
end