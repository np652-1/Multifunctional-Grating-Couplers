function [] = plotFieldonStructure_show_sim_settings_custom_colormap(hP,xP,yP,xSim,ySim,field,density,abs_patch_x,abs_patch_y,gc_patch_x,gc_patch_y,flux_Pts,title_string,dirname1,visibility,show_sim_settings,custom_colormap,climit)
% This function overlay fields on structure
% field: the quantity to be plotted, can be field amplitude or field
% intensity
% custom_colormap: the colormap used to plot the fields, use meep for field
% amplitude
% climit: 1 by 2 vector specifying the lower and upper limit for the color
% plot.

density_plot = interp2(xSim, ySim, density, xP, yP, 'nearest', 0); % 'nearest' is the best 
fig = figure('Visible',visibility); % select whether to show the image

ax1 = axes; % Plot the structure on the first axis
imagesc(ax1,xP(1,:), yP(:,1), max(density_plot,[],'all')-density_plot); % plot structure on axis 1
colormap(ax1, 'gray');
ax1.Visible='on'; % Make sure the first axis is visible
ax1.YDir = 'normal'; % set the direction of Y axis to be normal
% ax1.FontSize = 3; % set font size
daspect(ax1,[1,1,1])

% Plot the fields on the second, overlaid axis
ax2 = axes;
imagesc(ax2,xP(1,:), yP(:,1), field); % overlay fields on structure, here I plot the intensity
%c = meep;
colormap(ax2, custom_colormap);
colorbar(ax2,'Location','eastoutside')

% Adjust transparency of the second image
set(ax2.Children, 'AlphaData', 0.7); % 70% transparency for second image

% Set the second axis to be transparent and overlayed on the first
ax2.Color = 'none';
ax2.XTick = [];
ax2.YTick = [];
ax2.Position = ax1.Position; % Align the axes positions
set(ax2,'clim',climit)

ax2.UserData = linkprop([ax1,ax2],... % link the two axes
    {'Position','InnerPosition','DataAspectRatio', ...
    'ydir','xdir','xlim','ylim'}); % add more props as needed

% add a title
title(ax2,title_string)

if show_sim_settings=="on"
    % add patches to show the absorbers, gc region, and monitors
    hold(ax2,"on")
    x_min = min(min(xP)); % ref point to align patch with grids
    y_min = min(min(yP)); % ref point to align patch with grids
    patch(ax2, round((abs_patch_x-x_min)/hP)*hP+x_min-1/2*hP,round((abs_patch_y-y_min)/hP)*hP+y_min-1/2*hP,'w','FaceAlpha',0.1,'EdgeColor','k','linewidth',0.5) % show absorbers, I use round to match the grid, -1/2*hP is needed because the image pixels are centered on grid positions
    patch(ax2, round((gc_patch_x-x_min)/hP)*hP+x_min-1/2*hP,round((gc_patch_y-y_min)/hP)*hP+y_min-1/2*hP,'w','FaceAlpha',0.1,'EdgeColor','g','linewidth',0.5) % show gc region, I use round to match the grid, -1/2*hP is needed because the image pixels are centered on grid positions
    line(flux_Pts(:,1),flux_Pts(:,2),'color','red','linewidth',1); % show monitor lines
    hold(ax2,"off")
end

% export figure
exportgraphics(fig,[dirname1, '.jpg'],'Resolution',300); % use fig handle to show all axes, if you use gca, only ax2 would be exported.
end