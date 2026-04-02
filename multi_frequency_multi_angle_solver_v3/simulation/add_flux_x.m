function [flux_Pts,Ez_inc,Hy_inc] = add_flux_x(flux_center,flux_length, h, Ez_inc_fx,Hy_inc_fx,k,theta)
% Create a flux monitor along vertical line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% flux_center: center position of flux monitor, [x_pos,y_pos]
% flux_length: length of the flux monitor
% h: grid size on the flux monitor, I set it to be the same as simulation
% grid size
% Ez_inc: function handle that returns the z component of the incident 
% electric field on flux_Pts at a given k and theta
% Ez_inc: function handle that returns the y component of the incident 
% magnetic field on flux_Pts at a given k and theta
% k: wavenumber, 2*pi/wavelength, size (1, Nk)
% theta: angles of incident Gaussian beam, size (Nk, Nth)

% Outputs:
% flux_Pts: [x(:), y(:)] of points on flux monitor, size (N,2)
% Ez_inc:  z component of the incident electric field on flux monitor, size
% (N, Nk, Nth)
% Hy_inc:  y component of the incident magnetic field on flux monitor, size
% (N, Nk, Nth)
% Notes:
% monitor length can't be very large, otherwise we're overestimating the 
% power inside the waveguide by collecting the power of radiating modes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_pos = shift_if_on_grid(flux_center(1),h); % make sure x_pos is not on simulation grid to avoid singularities
y_pos = (flux_center(2)-flux_length/2):h:(flux_center(2)+flux_length/2); % y position of points on monitor, grid size h
flux_Pts = [repmat(x_pos,length(y_pos),1), y_pos(:)]; 
Nk = length(k);
Nth = size(theta,2); % Nth is the number of pitches, 1 wl 1 pitch determins 1 angle
Ez_inc = zeros(length(y_pos),Nk,Nth);
Hy_inc = zeros(length(y_pos),Nk,Nth);
for i=1:Nth % Nth is smaller than Nk, put it on the outer loop
    for j=1:Nk
        Ez_inc(:,j,i) = Ez_inc_fx(flux_Pts, k(j), theta(j,i)); % Incident electric field Ez on the flux monitor
        Hy_inc(:,j,i) = Hy_inc_fx(flux_Pts, k(j), theta(j,i)); % incident magnetic field Hy on the flux monitor
    end
end
end