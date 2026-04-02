function [CE, F1, F2, Einc_adj] = coupling_efficiency_fom(k, xySim, P, hs, Ezinc_fx, recPts, hr, P_wg, Ht_x, Pt)
% Calculate the coupling efficiency and adjoint incident field
% The integral is evaluated along a horizontal line above the waveguide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: wavenumber, 2*pi/wavelength, must be a scalar
% xySim: x y points in the simulation region, size (N1,2)
% P: polarization fields in the simulation region, size (ny, nx), no k^2
% prefactor, N1 = ny*nx
% hs: grid size in the source region
% Ezinc_fx: returns incident field Ez at a point (x, y)
% recPts: monitor points, size (N2,2), recPts(:,1) x, recPts(:,2) y
% hr: grid size on the monitor
% P_wg: power in the waveguide, scalar
% Ht_x: target magnetic field profile Hx, size (N2,1)
% Pt: power in the target mode

% Outputs:
% CE: coupling efficiency, ref: Zhexin Zhao, et al, Journal of Lightwave
% Technolohy, 2020
% F1: real part of the overlap integral
% F2: imaginary part of the overlap integral
% F1 and F2 are needed for gradients calculation
% Einc_adj: incident fields for adjoint simulation, size (ny, nx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ez_scat = electric_dipole_z_field_Ez(k, recPts, xySim, hs^2, P);
Ez_inc = Ezinc_fx(recPts(:,1),recPts(:,2));
Ez = Ez_scat + Ez_inc;
F1 = 1/2*real(sum(Ez.*conj(Ht_x)*hr))/sqrt(P_wg*Pt); % real part of the overlap integral
F2 = 1/2*imag(sum(Ez.*conj(Ht_x)*hr))/sqrt(P_wg*Pt); % imaginary part of the overlap integral
CE = F1^2+F2^2;

pz_adj = 1/4*conj(Ht_x)/sqrt(P_wg*Pt); % adjoint source dipole amplitude
Einc_adj = electric_dipole_z_field_Ez(k, xySim, recPts, hr, pz_adj); % Ez from pz_adj
Einc_adj = reshape(Einc_adj, size(P)); % adjoint incident field
end

