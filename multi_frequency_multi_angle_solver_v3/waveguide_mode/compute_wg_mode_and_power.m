function [Hy_t, P_t] = compute_wg_mode_and_power(k,h,mode_indices,flux_Pts,flux_center,flux_dir,wgt,n_core)
% Calcualtes target waveguide modes and their power
% flux_dir: +1 for positive and -1 for negative
% k: wave vector, size (Nk,1)
% mode_indices: waveguide mode indices, starts from 0, size (NK, Nth)
% At each wavelength, we calculate Nth modes. 
Nk = length(k);
Nth = size(mode_indices,2);
Hy_t = zeros(size(flux_Pts,1),Nk,Nth);
P_t = zeros(Nk,Nth);
for i=1:Nth
    for j=1:Nk
        [~,Ez_fx,Hy_fx] = TE_symmetric_waveguide(wgt,k(j),n_core,1,mode_indices(j,i)); % n_padding=1
        Hy_t(:,j,i) = Hy_fx(flux_Pts(:,2)-flux_center(2)); % These functions assume the waveguide's y_center is at 0
        Ez_t = Ez_fx(flux_Pts(:,2)-flux_center(2)); % These functions assume the waveguide's y_center is at 0
        P_t(j,i) = flux_dir*1/2*sum(real(Ez_t.*conj(Hy_t(:,j,i))*h)); % waveguide mode propagate along x, Px>0 
    end
end
end