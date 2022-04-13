% Mitchell Chandler, SIO
% Last updated: 12/04/2021

% >> run each of these after running the script for that transect

%% IX21
figure('color','w')
hold on
box on

%land mask:
contourf(topo_long2,topo_lat2,topo_bath2,[0 0],'LineColor','none')
colormap(0.7*[1 1 1])

%bathymetry:
deep_cmap = cmocean('deep');
contour(topo_long2,topo_lat2,topo_bath2,[-200 -200],'LineColor',deep_cmap(40,:),'LineWidth',1.5) %200 m isobath
contour(topo_long2,topo_lat2,topo_bath2,[-1000 -1000],'LineColor',deep_cmap(100,:),'LineWidth',1.5) %1000 m isobath
contour(topo_long2,topo_lat2,topo_bath2,[-2000 -2000],'LineColor',deep_cmap(200,:),'LineWidth',1.5) %2000 m isobath

%depth-integrated velocities:
scale = 0.02;
%N-ward vectors
pos = depth_int_v>=0; 
quiver(long_nom(pos),lat_nom(pos),scale*depth_int_u(pos),scale*depth_int_v(pos),'Color',[255 160 0]/255,'LineWidth',1,'Autoscale','off','ShowArrowHead','off')
%S-ward vectors
neg = depth_int_v<0; 
quiver(long_nom(neg),lat_nom(neg),scale*depth_int_u(neg),scale*depth_int_v(neg),'Color',[151 0 255]/255,'LineWidth',1,'Autoscale','off','ShowArrowHead','off')

xlim([30 60])
ylim([-36 -21])
daspect([1 1 1])

%axis ticks and labels:
XX = [30:5:60];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-36:5:-21];
yticks(YY)
YT = compose('%.0f\\circS',-YY);
yticklabels(YT)
set(gca,'TickDir','out','FontSize',12)

%% PX30
figure('color','w')
% clf
hold on
box on

%land mask:
contourf(topo_long2,topo_lat2,topo_bath2,[0 0],'LineColor','none')
colormap(0.7*[1 1 1])

%bathymetry:
deep_cmap = cmocean('deep');
contour(topo_long2,topo_lat2,topo_bath2,[-200 -200],'LineColor',deep_cmap(40,:),'LineWidth',1.5) %200 m isobath
contour(topo_long2,topo_lat2,topo_bath2,[-1000 -1000],'LineColor',deep_cmap(100,:),'LineWidth',1.5) %1000 m isobath
contour(topo_long2,topo_lat2,topo_bath2,[-2000 -2000],'LineColor',deep_cmap(200,:),'LineWidth',1.5) %2000 m isobath

%depth-integrated velocities:
scale = 0.02;
%N-ward vectors
pos = depth_int_v>=0; 
quiver(long_nom(pos),lat_nom(pos),scale*depth_int_u(pos),scale*depth_int_v(pos),'Color',[255 160 0]/255,'LineWidth',1,'Autoscale','off','ShowArrowHead','off')
%S-ward vectors
neg = depth_int_v<0; 
quiver(long_nom(neg),lat_nom(neg),scale*depth_int_u(neg),scale*depth_int_v(neg),'Color',[151 0 255]/255,'LineWidth',1,'Autoscale','off','ShowArrowHead','off')

xlim([150 180])
ylim([-33 -18])
daspect([1 1 1])

%axis ticks and labels:
XX = [150:5:180];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-33:5:-18];
yticks(YY)
YT = compose('%.0f\\circS',-YY);
yticklabels(YT)
set(gca,'TickDir','out','FontSize',12)

%% PX40
figure('color','w')
clf
hold on
box on

%land mask:
contourf(topo_long2,topo_lat2,topo_bath2,[0 0],'LineColor','none')
colormap(0.7*[1 1 1])

%bathymetry:
deep_cmap = cmocean('deep');
contour(topo_long2,topo_lat2,topo_bath2,[-200 -200],'LineColor',deep_cmap(40,:),'LineWidth',1.5) %200 m isobath
contour(topo_long2,topo_lat2,topo_bath2,[-1000 -1000],'LineColor',deep_cmap(100,:),'LineWidth',1.5) %1000 m isobath
contour(topo_long2,topo_lat2,topo_bath2,[-2000 -2000],'LineColor',deep_cmap(200,:),'LineWidth',1.5) %2000 m isobath

%depth-integrated velocities:
scale = 0.02;
%N-ward vectors
pos = depth_int_v>=0; 
quiver(long_nom(pos),lat_nom(pos),scale*depth_int_u(pos),scale*depth_int_v(pos),'Color',[255 160 0]/255,'LineWidth',1,'Autoscale','off','ShowArrowHead','off')
%S-ward vectors
neg = depth_int_v<0; 
quiver(long_nom(neg),lat_nom(neg),scale*depth_int_u(neg),scale*depth_int_v(neg),'Color',[151 0 255]/255,'LineWidth',1,'Autoscale','off','ShowArrowHead','off')
%Scale arrow
quiver(158.5,37,scale*0',scale*200','Color','k','LineWidth',2,'Autoscale','off','ShowArrowHead','off')
text(159,39.2,'200 m^2/s','Color','k','FontWeight','bold','FontSize',12)

xlim([135 165])
ylim([28 43])
daspect([1 1 1])

%axis ticks and labels:
XX = [135:5:165];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [28:5:43];
yticks(YY)
YT = compose('%.0f\\circN',YY);
yticklabels(YT)
set(gca,'TickDir','out','FontSize',12)

%label Izu Ridge
t = text(137.6,28.5,'Izu Ridge','Color','k','FontSize',12);
set(t,'Rotation',90);


