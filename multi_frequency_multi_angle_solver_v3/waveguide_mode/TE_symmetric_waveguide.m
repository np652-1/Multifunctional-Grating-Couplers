function [beta,Ez_fx,Hy_fx] = TE_symmetric_waveguide(wgt,k,n_core,n_padding,m)
% Calculates TE propagation constant and field profiles for a symmetric
% waveguide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% wgt: waveguide thickness
% k: vacuum wave number
% n_core: refractive index of waveguide core material
% n_padding: refractive index of waveguide padding material
% m: mode index, 0,1,2,3.....

% Outputs:
% beta: propagation constant
% Ez_fx: function that takes y coordinates and returns Ez fields
% Hy_fx: function that takes y coordinates and returns Hy fields

% Notes:
% The waveguide is along x direction, z is out-of-plane direction
% ref: Fundamentals of optical waveguides, Katsunari Okamoto, 2005, Ch2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = wgt/2; % half waveguide thickness
v = k*a*sqrt(n_core^2-n_padding^2); % normalized frequency
m_max = floor(v*2/pi); % maximum mode number at this wl, starts with 0
if m>m_max % mode m does not exist
    beta = NaN;
    Ez_fx = @(x) NaN;
    Hy_fx = @(x) NaN;
else
    TE_symmetric_dispersion = @(b) v*sqrt(1-b) - m*pi/2 - atan(sqrt(b/(1-b)));
    b = fzero(TE_symmetric_dispersion,[0,1]); % b belongs to [0,1]
    u = v*sqrt(1-b);
    w = v*sqrt(b);
    beta = sqrt(b*(n_core^2-n_padding^2)+n_padding^2)*k; % propagation constant
    kappa = u/a;
    ksi = w/a;
    phi = m*pi/2; % phase
    
    Ez_fx = @(x) (x>a).*cos(kappa*a-phi).*exp(-ksi*(x-a)) + (x<-a).*cos(kappa*a+phi).*exp(ksi*(x+a)) + ...
                    (x>=-a & x<=a).*cos(kappa*x-phi); % Ez field
    Hy_fx = @(x) -beta/k*((x>a).*cos(kappa*a-phi).*exp(-ksi*(x-a)) + (x<-a).*cos(kappa*a+phi).*exp(ksi*(x+a)) + ...
                    (x>=-a & x<=a).*cos(kappa*x-phi)); % Hy = -beta/k*Ez, mu0 = epsilon0 = 1
end