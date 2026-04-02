function [P_ext, Einc_adj] = power_extinction_fom(k, h, P, Einc, chi, Cpp, norm)
% Calculates power extinction and Einc for adjoint simulation
% Extinction is an explicit function of susceptibility, and the gradients
% have an extra term 1/2*k*h^2*imag(E.*conj(Einc))/norm besides 2*h^2*real(E_dir.*E_adj)
P_ext = 1/2*k*sum(sum(imag(P.*conj(Einc))))*h^2/norm;
pz_adj = -1j/4*k*chi.*conj(Einc)/norm;
einc_adj = BTTB_matvec_1(Cpp, k^2*pz_adj); 
Einc_adj = reshape(einc_adj, size(Einc));
end