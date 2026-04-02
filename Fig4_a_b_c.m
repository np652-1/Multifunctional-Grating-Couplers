%%
addpath('export_fig/');
addpath('tradeoff_curve_data/');

thickness_h_list = [10,12,16,20,24,32,36,40];
num_dp = 0;
plw=8;
fig=figure('units','normalized','position',[0,0,1,1]);
colors=winter(length(thickness_h_list));
hold on
for i2=1:length(thickness_h_list)
    thickness_h=thickness_h_list(i2);
    temp = readmatrix(fullfile(['tradeoff_curve_data/4wl_fom_opts_vs_num_layer_',num2str(thickness_h),'h.txt']));
    plot(temp(:,1),temp(:,2),'color',colors(i2,:),'LineWidth',plw)
    num_dp = num_dp + length(temp(:,1));
end

for i2=1:length(thickness_h_list)
    thickness_h=thickness_h_list(i2);
    temp = readmatrix(fullfile(['tradeoff_curve_data/4wl_fom_opts_vs_num_layer_',num2str(thickness_h),'h.txt']));
    scatter(temp(:,1),temp(:,2),500,colors(i2,:),'filled')
end
hold off
grid off
set(gca,'fontsize',35)

xticks([0:2:20])
yticks(0:0.2:1.0)
xlims=[0,20];
xlim(xlims)
ylims=[0.5,1];
ylim(ylims)
daspect([(xlims(2)-xlims(1)),2.5*(ylims(2)-ylims(1)),1])
export_fig('trade_off_4wls_nogrid.png','-transparent',fig)

%%
fig=figure('units','normalized','position',[0,0,1,1]);
hold on
for i2=1:8
num_layer=readmatrix(['tradeoff_curve_data/8wl_',num2str(i2),'_num_layer.txt']);
ce=readmatrix(['tradeoff_curve_data/8wl_',num2str(i2),'_fom_opts.txt']);
plot(num_layer,ce,'color',colors(i2,:),'LineWidth',plw)
end
for i2=1:8
num_layer=readmatrix(['tradeoff_curve_data/8wl_',num2str(i2),'_num_layer.txt']);
ce=readmatrix(['tradeoff_curve_data/8wl_',num2str(i2),'_fom_opts.txt']);
scatter(num_layer,ce,500,colors(i2,:),'filled')
end
%}
hold off
xticks([0:2:100])
yticks([])
grid off
xlims=[0,24];
xlim(xlims)
ylims=[0.5,1];
ylim(ylims)
set(gca,'fontsize',35)
daspect([(xlims(2)-xlims(1)),2.5*(ylims(2)-ylims(1)),1])
export_fig('trade_off_8wls_nogrid.png','-transparent',fig)