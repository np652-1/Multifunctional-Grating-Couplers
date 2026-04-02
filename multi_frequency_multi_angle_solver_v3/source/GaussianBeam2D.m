function [Einc] = GaussianBeam2D(Pts, k, theta, w0, z0, x0, y0)
% This function produces a 2D GaussianBeam as the incident electric field
% Ref: https://doc.comsol.com/6.1/doc/com.comsol.help.roptics/roptics_ug_optics.6.62.html
% k: wave vector
% w0: waist radius
% z0: axial distance from beam focus to center of illuminated area on the
% upper surface
% theta: angle between the propagation direction of the beam and -y axis 
% (-y ccw rotate to the propagation direction of the beam)
% the intersection of the beam's axis and the upper surface of the design
% region is (x0, y0)
% The fiels are evaluated at Pts
% reference: https://en.wikipedia.org/wiki/Gaussian_beam
x = Pts(:,1);
y = Pts(:,2);
z = z0 + ((x-x0)*tan(theta)+(y0-y))*cos(theta);
r = abs((x-x0-(y0-y)*tan(theta))*cos(theta));
zR = k*w0^2/2;
w = w0*sqrt(1+(z/zR).^2);
R = z.*(1+(zR./z).^2);
psi = 1/2*atan(z/zR); % add a 1/2 prefactor 
phase = exp(1j*(k*z+k*r.^2./(2*R)-psi));
phase(isnan(phase)) = 1; % when z=0 R=infinity phase = 1, get rid of singularities
Einc = sqrt(w0./w).*exp(-r.^2./w.^2).*phase; % E0 = 1, unit amplitude at focus, sqrt root of w0/w
end
