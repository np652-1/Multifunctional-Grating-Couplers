function [Ez_inc_fx,Hy_inc_fx,Einc,P_inc] = GaussianBeam2D_EzHy(k,theta,w0,z0,xySim,xSim,ySim,x,h,beam_center)
% Calcualate incident fields and power, provide Ez Hy function handles to 
% evaluate fields on given position 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: wavenumber, 2*pi/wavelength, size (1, Nk)
% theta: angles of incident Gaussian beam, size (Nk, Nth)
% w0: waist radius of Gaussian beam, scalar
% z0: distance between beam focus on the reference points on the structure
% xySim: [x(:),y(:)] for points in the simulation region, size (ny*nx,2)
% xSim, ySim: grid mesh, size (ny,nx)
% x: x cooridnates for points in the simulation region, size (1,nx)
% h: grid size
% beam_center: the intersection of the beam's axis and the structure's
% upper surface

% Outputs:
% Einc: incident electric fields Ez, size (ny,nx,Nk,Nth)
% P_inc: incident power, here I assumed the Gaussian beam is TEM, and this 
% approximation holds for wide beams (w0>>wavelength), size (Nk,Nth)
% Ez_inc: a function that returns Ez incident field on given points at
% certain k and theta
% Hy_inc: a function that returns Hy incident field on given points at
% certain k and theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ez_inc_fx = @(Pts,k, theta) GaussianBeam2D(Pts, k, theta, w0, z0, beam_center(1), beam_center(2)); 
Hy_inc_fx = @(Pts,k,theta) -1*GaussianBeam2D(Pts, k, theta, w0, z0, beam_center(1), beam_center(2))*sin(theta); % (ax/2, wg_y2) is the center of the top surface of the design region, Hy, Gaussian beam is TEM, and mu0=epsilon0=1 in the solver
Nk = length(k);
Nth = size(theta,2); % Nth is the number of pitches, 1 wl 1 pitch determins 1 angle
[ny,nx] = size(xSim);
Einc = zeros(ny,nx,Nk,Nth);
P_inc = zeros(Nk,Nth);
for i=1:Nth % Nth is smaller than Nk, put it on the outer loop
    for j=1:Nk
        Einc(:,:,j,i) = reshape(Ez_inc_fx(xySim,k(j),theta(j,i)), ny,nx); % incident field
        E0 = interp2(xSim, ySim, Einc(:,:,j,i), x, beam_center(2)); % get E on top surface
        P_inc(j,i) = sum(abs(E0).^2*h)/2*cos(theta(j,i)); % get incident power, this is the total power in the beam
        fprintf('\nIncident power at wavelength %g um and angle %g deg: %g \n', 2*pi/k(j),rad2deg(theta(j,i)),P_inc(j,i)); 
    end
end