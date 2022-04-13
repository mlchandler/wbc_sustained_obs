% Mitchell Chandler, SIO
% Last updated: 13/04/2022

%Have used colours from Paul Tol to ensure colourblind friendly palette https://personal.sron.nl/~pault/

load ix21_variability.mat
load px30_variability.mat
load px40_variability.mat

%% Plot annual cycle
figure('color','w')
clf

% --Agulhas--
%wbc transport
subaxis(3,3,1,'SpacingVert',2.5E-2)
hold on
shaded_error(1:12,ix21_wbc_transport_month,ix21_wbc_transport_month_SE,[0 0 0],0.1,4)
yline(mean(ix21_wbc_transport_raw),'--','LineWidth',2)
ylabel({'Western','boundary current','transport [Sv]'},'FontWeight','bold')
xlim([1 12])
xticks([1:12])
xticklabels([])
box on
title(['Agulhas (',num2str(abs(ix21_core_lat)),'^\circS)'],'FontWeight','bold')
set(gca,'FontSize',15)
ylim([-51 -39])
YL = ylim;
ypos = max(YL) - range(YL)*0.1;
text(1.25,ypos,'(a)','FontSize',15,'FontWeight','bold')
%core speed
subaxis(3,3,4)
hold on
shaded_error(1:12,ix21_core_speed_month,ix21_core_speed_month_SE,[0 68 136]/256,0.1,4)
yline(mean(ix21_core_speed_raw),'--','LineWidth',2)
ylabel({'Core','depth-integrated','velocity [m^2/s]'},'FontWeight','bold')
xlim([1 12])
xticks([1:12])
xticklabels([])
box on
set(gca,'FontSize',15)
ylim([-350 -290])
YL = ylim;
ypos = max(YL) - range(YL)*0.1;
text(1.25,ypos,'(d)','FontSize',15,'FontWeight','bold')
%offshore edge deviations
subaxis(3,3,7)
hold on
shaded_error(1:12,ix21_offshore_deviations_month,ix21_offshore_deviations_month_SE,[221 170 51]/256,0.2,4)
yline(mean(ix21_offshore_dev_raw),'--','LineWidth',2)
ylabel({'Offshore edge','deviations [km]'},'FontWeight','bold')
xlim([1 12])
xticks([1:12])
xticklabels(['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'])
box on
set(gca,'FontSize',15)
ylim([-40 40])
YL = ylim;
ypos = max(YL) - range(YL)*0.1;
text(1.25,ypos,'(g)','FontSize',15,'FontWeight','bold')

% --EAC--
%wbc transport
subaxis(3,3,2,'SpacingVert',2.5E-2)
hold on
shaded_error(1:12,px30_wbc_transport_month,px30_wbc_transport_month_SE,[0 0 0],0.1,4)
yline(mean(px30_wbc_transport_raw),'--','LineWidth',2)
xlim([1 12])
xticks([1:12])
xticklabels([])
box on
title(['EAC (',num2str(abs(px30_core_lat)),'^\circS)'],'FontWeight','bold')
set(gca,'FontSize',15)
ylim([-26 -16])
YL = ylim;
ypos = max(YL) - range(YL)*0.1;
text(1.25,ypos,'(b)','FontSize',15,'FontWeight','bold')
%core speed
subaxis(3,3,5)
hold on
shaded_error(1:12,px30_core_speed_month,px30_core_speed_month_SE,[0 68 136]/256,0.1,4)
yline(mean(px30_core_speed_raw),'--','LineWidth',2)
xlim([1 12])
xticks([1:12])
xticklabels([])
box on
set(gca,'FontSize',15)
ylim([-280 -215])
YL = ylim;
ypos = max(YL) - range(YL)*0.1;
text(1.25,ypos,'(e)','FontSize',15,'FontWeight','bold')
%offshore edge deviations
subaxis(3,3,8)
hold on
shaded_error(1:12,px30_offshore_deviations_month,px30_offshore_deviations_month_SE,[221 170 51]/256,0.2,4)
yline(mean(px30_offshore_dev_raw),'--','LineWidth',2)
xlim([1 12])
xticks([1:12])
xticklabels(['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'])
box on
set(gca,'FontSize',15)
ylim([-30 30])
YL = ylim;
ypos = max(YL) - range(YL)*0.1;
text(1.25,ypos,'(h)','FontSize',15,'FontWeight','bold')

% --Kuroshio--
%wbc transport
subaxis(3,3,3,'SpacingVert',2.5E-2)
hold on
shaded_error(1:12,px40_wbc_transport_month,px40_wbc_transport_month_SE,[0 0 0],0.1,4)
yline(mean(px40_wbc_transport_raw),'--','LineWidth',2)
xlim([1 12])
xticks([1:12])
xticklabels([])
box on
title(['Kuroshio (',num2str(abs(px40_core_lat)),'^\circN)'],'FontWeight','bold')
set(gca,'FontSize',15)
ylim([61 79])
YL = ylim;
ypos = max(YL) - range(YL)*0.1;
text(1.25,ypos,'(c)','FontSize',15,'FontWeight','bold')
%core speed
subaxis(3,3,6)
hold on
shaded_error(1:12,px40_core_speed_month,px40_core_speed_month_SE,[0 68 136]/256,0.1,4)
yline(mean(px40_core_speed_raw),'--','LineWidth',2)
xlim([1 12])
xticks([1:12])
xticklabels([])
box on
set(gca,'FontSize',15)
ylim([400 540])
YL = ylim;
ypos = max(YL) - range(YL)*0.1;
text(1.25,ypos,'(f)','FontSize',15,'FontWeight','bold')
%offshore edge deviations
subaxis(3,3,9)
hold on
shaded_error(1:12,px40_offshore_deviations_month,px40_offshore_deviations_month_SE,[221 170 51]/256,0.2,4)
yline(mean(px40_offshore_dev_raw),'--','LineWidth',2)
xlim([1 12])
xticks([1:12])
xticklabels(['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'])
box on
set(gca,'FontSize',15)
ylim([-40 40])
YL = ylim;
ypos = max(YL) - range(YL)*0.1;
text(1.25,ypos,'(i)','FontSize',15,'FontWeight','bold')


