function [CE, F1, F2, Einc_adj] = coupling_efficiency_fom_v2(k, xySim, P, hs, Ez_inc, recPts, hr, P_inc, Hy_t, Pt)
% Calculate the coupling efficiency and adjoint incident field
% The integral is evaluated along a line perpendicualr to the waveguide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: wavenumber, 2*pi/wavelength, must be a scalar
% xySim: x y points in the simulation region, size (N1,2)
% P: polarization fields in the simulation region, size (ny, nx), no k^2
% prefactor, N1 = ny*nx
% hs: grid size in the source region
% Ezinc: incident field Ez on monitor, size (N2,1)
% recPts: monitor points, size (N2,2), recPts(:,1) x, recPts(:,2) y
% hr: grid size on the monitor
% P_inc: incident power from the source, scalar
% Hy_t: target magnetic field profile Hy, size (N2,1)
% Pt: power in the target mode

% Outputs:
% CE: coupling efficiency, ref: Zhexin Zhao, et al, Journal of Lightwave
% Technolohy, 2020
% F1: real part of the overlap integral
% F2: imaginary part of the overlap integral
% F1 and F2 are needed for gradients calculation
% Einc_adj: incident fields for adjoint simulation, size (ny, nx)

% Notes:
% flux_dir does not matter for this fom because of the abs operation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ez_scat = electric_dipole_z_field_Ez(k, recPts, xySim, hs^2, P);
Ez = Ez_scat + Ez_inc;
F1 = 1/2*real(sum(Ez.*conj(Hy_t)*hr))/sqrt(P_inc*Pt); % real part of the overlap integral
F2 = 1/2*imag(sum(Ez.*conj(Hy_t)*hr))/sqrt(P_inc*Pt); % imaginary part of the overlap integral
CE = F1^2+F2^2;

pz_adj = 1/4*conj(Hy_t)/sqrt(P_inc*Pt); % adjoint source dipole amplitude
Einc_adj = electric_dipole_z_field_Ez(k, xySim, recPts, hr, pz_adj); % Ez from pz_adj
Einc_adj = reshape(Einc_adj, size(P)); % adjoint incident field
end

