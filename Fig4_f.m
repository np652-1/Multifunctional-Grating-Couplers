%%
clear;
addpath(genpath('export_fig/'))
addpath(genpath('multi_frequency_multi_angle_solver_v3/'))
mkdir('Fig4_f_plots/')

disp('Solving and plotting started. Please wait ~30 seconds for plots to start generating.')

% basic parameters
wl=linspace(0.536,0.609,2);
wl_max = max(wl);
h=1/64;
n_core=2.0;
wgt=64*h;
gcl=floor(10.7/h)*h;
par_flag = 0;
k=2*pi./wl;

% solver settings
tol = 1e-5; % used in bicgstab
maxit = 10; % used in bicgstab
digits(16); % precision settings for matlab
quad_order = 6; % quadrature correction order, deals with Greens function sigularities

% simulation region
absl = 5*wl_max; % absorber length
padding = 3*wl_max; % distance between flux monitor and gc region.
wgl = gcl+2*absl+padding+2*wl_max; % waveguide length, include absorbers on two sides and paddings betwenn gc and flux monitor, and 1*wl_max between flux monitor and left absorber, 1*wl_max between gc region and right absorber. 
[ax,ay,Npv,Nph,nx,ny,x,y,xSim,ySim,xySim,NODES,param] = build_simulation(k,wgl,wgt,h,par_flag);

% waveguide region
wg_center = [ax/2,ay/2];
wg_size = [wgl,wgt];
[wg_x1,wg_x2,wg_y1,wg_y2,wg_region_index] = add_region(h,wg_center,wg_size,xSim,ySim);
%fprintf('\nWaveguide thickness: %g\n', wg_y2-wg_y1);

% set material and initialize chi
density = zeros(size(xSim)); % We don't need to add smoothed backgrounds for this, the error is very small
density(wg_region_index) = 1; % initial density matrix
chi_max = n_core^2-1; % permittivity
chi = chi_max*density; % initial susceptibility matrix

% Add absorbers on two sides
target_R = 1e-20; % Smaller value for smaller reflection
[chi,abs_patch_x,abs_patch_y] = add_absorber(2*pi/wl_max,target_R,absl,n_core,xSim,ySim,wg_x1,wg_x2,wg_y1,wg_y2,h,chi); % set absorbers according to min(k)

% grating coupler region, put it at the center
etch = wg_y2-wg_y1; % cut through the waveguide, wg_y2-wg_y1 is the actual size
gc_center = [wg_x2-absl-1*wl_max-gcl/2,wg_y2-etch/2]; % put the gc on the right, 1 wl_max away from the absorber
gc_size = [gcl,etch];
[gc_x1,gc_x2,gc_y1,gc_y2,gc_region_index,nx_gc,ny_gc,gc_patch_x,gc_patch_y] = add_region(h,gc_center,gc_size,xSim,ySim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Plot Optimal Design Density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numlayers=4;
numreps=ny_gc/numlayers;
x_opt = readmatrix('tradeoff_curve_data/b_gc_density_opt_betaInf.txt');
xoptre=reshape(x_opt,ny_gc/numreps,nx_gc);

xoptreactual=zeros(ny_gc,nx_gc);
ind=1;
for i2=1:size(xoptre,1)
xoptreactual(ind:(ind+numreps-1),:)=repmat(xoptre(i2,:),numreps,1);
ind=ind+numreps;
end
xoptactual=xoptreactual(:);

fig=figure('units','normalized','position',[0,0,1,1],'visible','off');
imagesc(1-xoptreactual)
colormap('gray')
daspect([1,1,1])
set(gca,'ydir','normal')
xticks([])
yticks([])
export_fig('Fig4_f_plots/optimal_design_density.png','-transparent',fig)

gc_density = reshape(xoptreactual,ny_gc,nx_gc);
chi(gc_region_index) = chi_max*xoptactual;
density(gc_region_index)=xoptactual;
x_gc=reshape(xSim(gc_region_index),ny_gc,nx_gc);
y_gc=reshape(ySim(gc_region_index),ny_gc,nx_gc);

%%
z0 = 0;
marker_pitch=[1.6,2.2];
w0=3.0;
theta = -asin(2*pi./reshape(marker_pitch,1,[])./k(:));
beam_center = [gc_center(1),wg_y2];
[Ez_inc_fx,Hy_inc_fx,Einc,P_inc] = GaussianBeam2D_EzHy(k,theta,w0,z0,xySim,xSim,ySim,x,h,beam_center);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Plot Fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorwls=[450,470;620,650];
folder='Fig4_f_plots';
for kind=1:2
    for thetaind=1:2
[afun, Chi, precon, Cpp] = build_solver(k(kind), chi, ax, ay, nx, ny, Npv, Nph, NODES{kind}, param{kind}, quad_order); % all angles at this freq share the outputs. NODES{1} is the NODES at k (single freq)
P = solve_P_E(k(kind), afun, Chi, precon, Cpp, Einc(:,:,kind,thetaind), tol, maxit); % forward simulation

LvP = ay+10*wl_max;
LhP = ax+5*wl_max;
hP = h; % hP can be different from h
ctrPt = [ax/2, ay/2+2];
Einc_fxP = @(x, y) Ez_inc_fx([x,y],k(kind),theta(kind,thetaind));
[EincP, EscatP, EP, xP, yP, title_string] = plotELarge(LvP, LhP, hP, k(kind), ctrPt, quad_order, xSim, ySim, P, Einc_fxP, ...
    ['\lambda=', num2str(wl(kind),'%.3f'),' um', ' \theta=', num2str(rad2deg(theta(kind,thetaind)),'%.1f'),' deg'],'off');

fig=figure('units','normalized','position',[0,0,1,1],'visible','off');
ax1=axes;%subplot('Position',[0,0,1,1]);
density_plot = interp2(xSim, ySim, density, xP, yP, 'nearest', 0); % 'nearest' is the best
imagesc(ax1,xP(1,:), yP(:,1), max(density_plot,[],'all')-density_plot); % plot structure on axis 1
colormap(ax1, 'gray');
ax1.Visible='on'; % Make sure the first axis is visible
ax1.YDir = 'normal'; % set the direction of Y axis to be normal
daspect(ax1,[1,1,1])
axis off
box off
xticks([])
yticks([])
xlims=[min(xP(1,:))+5*max(wl),max(xP(1,:))-6*max(wl)];
xlim(xlims)
ylims=[min(yP(:,1))+1.5*max(wl),max(yP(:,1))-5*max(wl)];
ylim(ylims)

% Plot the fields on the second, overlaid axis
ax2 = axes;
imagesc(ax2,xP(1,:), yP(:,1), abs(EP).^2); % overlay fields on structure, here I plot the intensity
colormap(ax2, singleColorIntensityColormap(colorwls(kind,thetaind)));
%colorbar(ax2,'Location','eastoutside')

% Adjust transparency of the second image
set(ax2.Children, 'AlphaData',0.85); % 70% transparency for second image
axis off
box off
ax2.Color = 'none';
ax2.XTick = [];
ax2.YTick = [];
ax2.Position = ax1.Position; % Align the axes positions
set(ax2,'clim',[0,max(abs(EP),[],'all')])
ax2.UserData = linkprop([ax1,ax2],... % link the two axes
    {'Position','InnerPosition','DataAspectRatio', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed
set(gca,'FontSize',22)

export_fig([folder,'/wl',num2str(kind),'_k',num2str(thetaind),'.png'],fig);
    end
end