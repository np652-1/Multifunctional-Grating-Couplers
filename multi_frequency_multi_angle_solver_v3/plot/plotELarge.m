function [EincP, EscatP, EP, xP, yP, title_string] = plotELarge(Lv, Lh, h, k, ctrPt, quad_order, xSim, ySim, P, Einc_fx, title_string,show_plot)
% Plot Einc, Escat, and E in a region larger than the simulation region
[ayP,~,nyP,~,axP,~,nxP,~] = identical_grid_info(Lv,Lh,h); % generate uniform grid 
[D0, D1]= planck_win_get_D(k,h);
CppP = build_greens_function(k, axP, ayP, nxP, nyP, quad_order, D0, D1); % get green's function

offset = [axP/2-ctrPt(1), ayP/2-ctrPt(2)]; % shift the plot region
[xP, yP] = meshgrid((0:h:axP)-offset(1), (0:h:ayP)-offset(2));
P_plot = interp2(xSim, ySim, P, xP, yP, 'linear', 0); % P is polarization strength in the simulation region
escatP = BTTB_matvec(CppP*k^2, P_plot);
EscatP = reshape(escatP, size(xP));
EincP = Einc_fx(xP(:),yP(:));
EincP =  reshape(EincP, size(xP));
EP = EincP + EscatP;

if show_plot=="on"
    c = meep;
    figure;
    imagesc(xP(1,:), yP(:,1), real(EincP));
    colormap(c);
    set(gca,'YDir','Normal')
    title([title_string ': Re(Einc)'])
    colorbar
    daspect([1 1 1])
    
    figure;
    imagesc(xP(1,:), yP(:,1), real(EscatP));
    colormap(c);
    set(gca,'YDir','Normal')
    title([title_string ': Re(Escat)'])
    colorbar
    daspect([1 1 1])
    
    figure;
    imagesc(xP(1,:), yP(:,1), real(EP));
    colormap(c);
    set(gca,'YDir','Normal')
    title([title_string ': Re(E)'])
    colorbar
    daspect([1 1 1])
end
end