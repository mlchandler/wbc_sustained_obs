% Mitchell Chandler, SIO
% Last updated: 08/10/2021

load ix21_fig2.mat
load px30_fig2.mat
load px40_fig2.mat

%% Plot
plot_range = 2:70;
yticklength = 2.5E-2;
xticklength = 8E-3;
fsize = 17;

mv = 600; %saturate colourbar

alt_optn = cmocean('balance');
alt_optn(127:130,:) = [1 1 1;1 1 1;1 1 1;1 1 1];

figure('color','w') 
clf

%IX21
subaxis(1,3,1,'SpacingHoriz',3E-2)
hold on
imagesc(ix21_long_nom(plot_range),time_monthly,ix21_divel(plot_range,:)')
contour(ix21_long_nom(plot_range),time_monthly,ix21_divel(plot_range,:)',[0 0],'Linecolor','w','LineWidth',1)
plot(ix21_long_nom(plot_range(end)),ix21_xbt_times,'wo','MarkerSize',10,'MarkerFaceColor','k') %timestamps for XBT transects
plot(ix21_core_long,time_monthly,'-','Color','y','LineWidth',3.5) %location of jet core
plot(ix21_offshore_long,time_monthly,'-','Color','k','LineWidth',3) %offshore edge
colormap(alt_optn)
caxis([-mv mv])
set(gca,'TickDir','out') 
ax=gca;
ax.YAxis.TickLength = [yticklength yticklength];
ax.XAxis.TickLength = [xticklength xticklength];
yticks(datenum(2005,1:24:204,1))
datetick('y','keeplimits','keepticks')
ylabel('Year')
xticks([ceil(ix21_long_nom(plot_range(1))):3:ix21_long_nom(plot_range(end))])
set(gca,'FontSize',fsize)
box off
title('(a) IX21, Agulhas Current')
  
%PX30
subaxis(1,3,2)
hold on
imagesc(px30_long_nom(plot_range),time_monthly,px30_divel(plot_range,:)')
contour(px30_long_nom(plot_range),time_monthly,px30_divel(plot_range,:)',[0 0],'Linecolor','w','LineWidth',1)
plot(px30_long_nom(plot_range(end)),px30_xbt_times,'wo','MarkerSize',10,'MarkerFaceColor','k') %timestamps for XBT transects
plot(px30_core_long,time_monthly,'-','Color','y','LineWidth',3.5) %location of jet core
plot(px30_offshore_long,time_monthly,'-','Color','k','LineWidth',3) %offshore edge
colormap(alt_optn)
caxis([-mv mv])
set(gca,'TickDir','out') 
ax=gca;
ax.YAxis.TickLength = [yticklength yticklength];
ax.XAxis.TickLength = [xticklength xticklength];
yticks(datenum(2005,1:24:204,1))
datetick('y','keeplimits','keepticks')
yticklabels([])
xticks([ceil(px30_long_nom(plot_range(1))):3:px30_long_nom(plot_range(end))])
xlabel('Longitude [\circE]')
set(gca,'FontSize',fsize)
box off
title('(b) PX30, EAC')

%PX40
subaxis(1,3,3)
hold on
imagesc(px40_long_nom(plot_range),time_monthly,px40_divel(plot_range,:)')
contour(px40_long_nom(plot_range),time_monthly,px40_divel(plot_range,:)',[0 0],'Linecolor','w','LineWidth',1)
plot(px40_long_nom(plot_range(end)),px40_xbt_times,'wo','MarkerSize',10,'MarkerFaceColor','k') %timestamps for XBT transects
plot(px40_core_long,time_monthly,'-','Color','y','LineWidth',3.5) %location of jet core
plot(px40_offshore_long,time_monthly,'-','Color','k','LineWidth',3) %offshore edge
colormap(alt_optn)
caxis([-mv mv])
set(gca,'TickDir','out') 
ax=gca;
ax.YAxis.TickLength = [yticklength yticklength];
ax.XAxis.TickLength = [xticklength xticklength];
yticks(datenum(2005,1:24:204,1))
datetick('y','keeplimits','keepticks')
yticklabels([])
xticks([ceil(px40_long_nom(plot_range(1))):3:px40_long_nom(plot_range(end))])
set(gca,'FontSize',fsize)
box off
title('(c) PX40, Kuroshio')

%Colorbar
c=colorbar;
set(c,'Position', [0.92, 0.2, 0.015, 0.6]);
ylabel(c,'Depth-integrated velocity [m^2/s]')
cbarrow






