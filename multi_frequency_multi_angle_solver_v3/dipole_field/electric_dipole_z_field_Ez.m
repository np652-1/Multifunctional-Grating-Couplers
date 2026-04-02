function [Ez] = electric_dipole_z_field_Ez(k, recPts, srcPts, srcdA, P)
% Calculate the electric field Ez from electric dipole Pz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% P: dipole amplitude Pz, size (N1,1) or (ny, nx) N1=ny*nx, no k^2 prefactor
% srcPts: source points, size (N1,2), srcPts(:,1) x, srcPts(:,2) y
% srcdA: source area, 2D source h^2, 1D source h, point source 1, h is grid size 
% recPts: field points, size (N2,2), recPts(:,1) x, recPts(:,2) y
% When recPts overlap with srcPts, G has singularities. When monitor is in
% the structure, we can put the monitor off grid so that distance is never
% zero.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = recPts(:,1); % field points
y = recPts(:,2);

x0 = srcPts(:,1); % source points
y0 = srcPts(:,2);

r = sqrt((x-x0.').^2 + (y-y0.').^2); % distance matrix

G = 1j/4*k^2*besselh(0,k*r)*srcdA;% Green's function matrix
% G(isnan(G)) = 0; % set the singularities to zero.
Ez = G*P(:);
end