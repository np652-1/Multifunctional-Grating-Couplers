function [NODES, param] = build_nodes_parfor(k, Lv, Lh, h, Npv, Nph, proxy_type)
% use parfor to speed up when there are a lot of frequencies.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% k: wavenumber, 2*pi/wavelength, size (N, 1)
% Lv: vertical length of the simulation region
% Lh: horizontal length of the simulation region
% h: grid size
% Npv: total vertical panel number, see identical_grid_info.m
% Nph: total horizontal panel number, see identical_grid_info.m
% proxy_type: select type of proxy layer, this affects calculation speed 

% Outputs:
% NODES: simulation nodes, size {N,1}, type cell, each element is a cell,
% which stores structured simulation info
% see identical_tree_generation for detailed data structure
% param: solver parameters and quadrature correction coefficients, size
% {N,1}, type cell, each element is a struct
    % param.quad_order: quadature order for direct solvers
    % Current code only support 2nd or 4th quadrature, choose eather 2 or 4
    % param.compress_acc: relative compression accuracy
    % param.proxy_layer: number of proxy layer, use 1 for compression acc > 1e-5
    % param.print_flag: print rank information during compression, 1 yes, 0 no
    % param.plot_flag: show compressed points plot for each level, 1 yes, 0 no
    % param.fileID: the place printouts go, 1 for the screen
    % param.D0: 4th order correction coefficients
    % param.D1: 6th order correction coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% direct solver parameters
param.quad_order = 4; % choose the quadature order for direct solvers
%%% Current code only support 2nd or 4th quadrature, choose eather 2 or 4
param.compress_acc = 1e-5; % relative compression accuracy
param.proxy_layer = 2; % number of proxy layer, use 1 for compression acc > 1e-5
param.print_flag = 0; % print rank information during compression, 1 yes, 0 no
param.plot_flag = 0; % show compressed points plot for each level, 1 yes, 0 no
param.fileID = 1; % the place printouts go, 1 for the screen

% create tree structure
x = 0:h:Lh;
y = 0:h:Lv;
NODES = identical_tree_generation(y,x,Npv,Nph);
% nlevels = NODES{1,end}; % number of levels of the binary tree

% Compute proxy surface for compression
switch proxy_type
    case 1
        NODES = near_surround_proxy_surface_creation(NODES,param.proxy_layer,Nph);
    case 2
        NODES = nonsurround_proxy_surface_creation(NODES,param.proxy_layer,Nph);
    case 3
        NODES = proxy_surface_creation(NODES,param.proxy_layer);
end
% store NODES and param for different k
if length(k)>1
    NODES = repmat({NODES}, length(k), 1);
    param = repmat({param}, length(k), 1);
    % myPool = parpool(feature('NumCores'));
    % myPool = parpool(8);
    tic
    parfor i = 1:length(k)
        [D0, D1]= planck_win_get_D(k(i),h); % get quadrature correction
        param{i}.D0 = D0; % store correction
        param{i}.D1 = D1;
        NODES{i} = identical_local_compression_all_proxy(NODES{i},k(i),h,param{i}); % compression
    end    
    tbuild = toc
    % delete(myPool);
else
    NODES = {NODES};
    param = {param};
    [D0, D1]= planck_win_get_D(k,h); % get quadrature correction
    param{1}.D0 = D0; % store correction
    param{1}.D1 = D1;
    tic
    NODES{1} = identical_local_compression_all_proxy(NODES{1},k,h,param{1}); % compression
    tbuild = toc
end
end