function [ax,ay,Npv,Nph,nx,ny,x,y,xSim,ySim,xySim,NODES,param] = build_simulation(k,lx,ly,h,par_flag)
% build simulation grid and NODES, NODES is used to calculate the inverse
% of (I-BG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: wavenumber, 2*pi/wavelength, size (N, 1)
% lx: horizontal size of the simulated structure
% ly: vertical size of the simulated structure
% h: grid size
% par_flag: set to 1 to use parfor

% Outputs:
% ax: horizontal size of simulation region (ax>lx)
% ay: vertical size of simulation region (ay>ly)
% Npv: total vertical panel number, see identical_grid_info.m
% Nph: total horizontal panel number, see identical_grid_info.m
% nx: number of grid points in horizontal direction
% ny: number of grid points in vertical direction
% x,y,xSim,ySim: grid points

% Notes:
% Use parfor to speed up when N is large
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Lh = lx+2*h; % make sure the simulation region is larger than the structure region
Lv = ly+2*h;
[ay,~,ny,Npv,ax,~,nx,Nph] = identical_grid_info(Lv,Lh,h); % generate uniform grid, the new domain is never smaller than the original
x = 0:h:ax; 
y = 0:h:ay; 
[xSim,ySim] = meshgrid(x,y);
xySim = [xSim(:) ySim(:)];
if par_flag==1
    [NODES, param] = build_nodes_parfor(k, ay, ax, h, Npv, Nph,1); % build simulation nodes, use NODES{index}, param{index} to access the data, we can use parfor to speed up
else
    [NODES, param] = build_nodes(k, ay, ax, h, Npv, Nph,1); % build simulation nodes, use NODES{index}, param{index} to access the data, we can use parfor to speed up
end
end