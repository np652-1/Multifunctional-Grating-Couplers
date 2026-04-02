function [CE, gradients] = compute_CE_and_gradients(density_design,gc_region_index,chi,chi_max,flux_Pts,Ezinc_monitor,Hy_t,Pt,Einc, P_inc,k,ax,ay,nx,ny,xySim,h,Npv,Nph,NODES,param,quad_order,tol,maxit)
% This function calculates coupling efficiency (overlap integral) and the gradients
% for all wavelengths and angles
chi(gc_region_index) = chi_max*density_design; % update the design region
Nk = size(P_inc,1);
Nth = size(P_inc,2);
CE = zeros(Nk,Nth);
gradients = zeros(length(density_design),Nk,Nth);
for i=1:Nk % use parfor to speed up
    [afun, Chi, precon, Cpp] = build_solver(k(i), chi, ax, ay, nx, ny, Npv, Nph, NODES{i}, param{i}, quad_order); % all angles at this freq share the outputs. NODES{i} is the NODES at k(i) 
    for j=1:Nth
        [P, E_dir] = solve_P_E(k(i), afun, Chi, precon, Cpp, Einc(:,:,i,j), tol, maxit); % forward simulation
        [CE(i,j), F1, F2, Einc_adj] = coupling_efficiency_fom_v2(k(i), xySim, P, h, Ezinc_monitor(:,i,j), flux_Pts, h, P_inc(i,j), Hy_t(:,i,j), Pt(i,j)); % calcualte CE and adjoint incident fields
        if nargout>1 % require gradients
            [~, E_adj] = solve_P_E(k(i), afun, Chi, precon, Cpp, Einc_adj, tol, maxit); % adjoint simulation
            gradients_sim = 2*F1*2*h^2*real(E_dir.*E_adj)*chi_max + 2*F2*2*h^2*real(E_dir.*(-1i*E_adj))*chi_max; % gradients in the simulation region
            gradients(:,i,j) = gradients_sim(gc_region_index); % gradients in the gc region
        end
        fprintf('\nCE at wavelength %g angle %g: %g\n\n',i,j,CE(i,j));
    end 
end