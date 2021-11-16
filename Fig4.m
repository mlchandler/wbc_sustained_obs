% Mitchell Chandler, SIO
% Last updated: 02/11/2021

%Colours from Paul Tol (https://personal.sron.nl/~pault/) and Brewermap to ensure colourblind friendly palettes 

load ix21_fig2.mat
load ix21_variability.mat
load ix21_velocity.mat

%% Find SLA points at Agulhas core along IX21
%Read in monthly SLA
sla_monthly = ncread('sla_monthly_2004-2019_indo-pac.nc','sla'); %referenced to 1993-2012 period [m]
lat_sat_monthly = ncread('sla_monthly_2004-2019_indo-pac.nc','latitude');
long_sat_monthly = ncread('sla_monthly_2004-2019_indo-pac.nc','longitude');

%Find mean long position of Agulas core
ix21_mean_core_long = mean(ix21_core_long);

%Find SLA lat/long grid points nearest to Agulhas mean core lat/long
[~,I_long] = min(abs(long_sat_monthly - ix21_mean_core_long));
long_band = long_sat_monthly(I_long);
[~,I_lat] = min(abs(lat_sat_monthly - ix21_core_lat));
lat_band = lat_sat_monthly(I_lat);
lat_N = lat_sat_monthly(I_lat+1);
lat_S = lat_sat_monthly(I_lat-1);

%Find SLA one grid point above (N) and below (S) the core for computing velocity 
sla_N = squeeze(sla_monthly(I_long,I_lat+1,:));
sla_S = squeeze(sla_monthly(I_long,I_lat-1,:));
%re-reference SLA to be relative to 2004-2019 mean
sla_N = sla_N - mean(sla_N);
sla_S = sla_S - mean(sla_S);

%% Compute geostrophic velocity
g=9.81;
f=gsw_f(lat_band);
deta = sla_N - sla_S;
dy = gsw_distance([long_band long_band],[lat_N lat_S]); %[m]
uvel = -g/f*deta/dy;

%% Correlate u' and IX21 core longitude
[r,p,eDOF] = corr_pval(uvel,ix21_core_long);
ix21_rpd = [r;p;eDOF]

%% Other correlations
% %Correlate u' and IX21 core depth-integrated velocity:
% [r,p] = corr_pval(uvel,ix21_core_speed_raw')
% 
% %Correlate u' and IX21 WBC transport:
% [r,p] = corr_pval(uvel,ix21_wbc_transport_raw)
% 
% %Correlate IX21 core longitude and core depth-integrated velocity:
% [r,p] = corr_pval(ix21_core_long,ix21_core_speed_raw')

%% Find 1 stdev of offshore deviations 
%Rouault and Penven 2011; Krug and Tournadre 2012; Elipot and Beal 2015
offshore_std = mean(ix21_core_long)+std(ix21_core_long);

%Index when offshore deviations exceed this
NP_idx = find(ix21_core_long >= offshore_std);

%Find number of Natal Pulse events
num_NP = numel(NP_idx);
av_num_NP = num_NP/16;
NP_events = [num_NP, av_num_NP]

%Re-reference SLA to be relative to 2004-2019 mean (for plotting)
sla_monthly0419 = sla_monthly - mean(sla_monthly,3);

%% Velocity cross-section composites for Natal Pulse and non-Natal Pulse periods
idx_long = find(ix21_long_nom < 37);

no_idx = find(ix21_core_long < offshore_std); %non-NP times

%composites
gvel_NP = mean(ix21_gvel_LKM(:,idx_long,NP_idx),3);
gvel_normal = mean(ix21_gvel_LKM(:,idx_long,no_idx),3);
gvel_diff = gvel_NP - gvel_normal;

%% Use bootstrapping to test whether the mean composites are significantly different
%This method treats each variability state as a separate population 
%and bootstraps from each population to produce samples A and B then
%computes the difference between them. The mean and CI for the bootstrapped
%differences is computed. The difference between the actual composites is
%determined to be significant if it falls within the bootstrapped CI, and
%if the boostrapped CI does not cross 0.

iterations = 1000; %1000 | 1E5
alpha = 0.05 %e.g. 0.1 = 90% CI, 0.05 = 95% CI

store_diff = NaN(size(gvel_diff,1),size(gvel_diff,2),iterations); %[depth x long x iterations]

for i=1:iterations
    %sample with replacement from NP state to build subsample of same size
    A = datasample(NP_idx,length(NP_idx));
    mean_A = mean(ix21_gvel_LKM(:,idx_long,A),3);
    
    %sample with replacement from non-NP state to build subsample of same size
    B = datasample(no_idx,length(no_idx));
    mean_B = mean(ix21_gvel_LKM(:,idx_long,B),3);
    
    %difference
    store_diff(:,:,i) = mean_A - mean_B;
end
  
%find upper and lower percentile values to give CI
lower_bound = prctile(store_diff,alpha/2,3);
upper_bound = prctile(store_diff,100-alpha/2,3);

%Composite difference is significant if within the 95% CI, and the 95% CI does not cross 0.
CI_mask = double(lower_bound.*upper_bound < 0 | gvel_diff < lower_bound | gvel_diff > upper_bound); %(not significant = 1)
% if the CI lower bound and upper bound have different signs then multiplying them will give a negative number
CI_mask(isnan(gvel_diff)) = NaN; %mask NaNs

%% -- Plot --
fsize = 13;
cticklength = 0.05;

figure('color','w')
clf

%Time series
subaxis(3,2,[1 2],'SpacingVert',6E-2,'SpacingHoriz',9E-2)
plot(time_monthly,ix21_core_long,'Color','k','LineWidth',3)
yline(offshore_std,'--','LineWidth',2)
ylim([32 34.5])
ylabel('Core longitude [\circE]')
xticks(datenum(2004,1:12*2:204,1))
xlim([datenum('01-Jan-2004') datenum('01-Jan-2020')])
datetick('x','keepticks','keeplimits')
box on
text(datenum('01-May-2004'),34.25,'(a)','FontWeight','bold','FontSize',fsize)
set(gca,'FontSize',fsize)
%yearly gridlines
D = datenum(2004,1:12:204,1);
for i=1:length(D)
xline(D(i),'Color',[0 0 0]+0.8)    
end


mv = 0.8;
%NP cross-section composite
subaxis(3,2,3)
hold on
set(gca,'Color','k') %background colour (and therefore NaNs)
contourf(ix21_long_nom(idx_long),-argo_depth,gvel_NP,-mv:0.1:mv,'LineColor','none') 
contour(ix21_long_nom(idx_long),-argo_depth,gvel_NP,[0 0],'LineColor',[0 0 0],'LineWidth',1.5) %0
contour(ix21_long_nom(idx_long),-argo_depth,gvel_NP,-mv:0.1:-0.1,'--','LineColor',[0 0 0]+0.5,'LineWidth',1.5) %-ve
contour(ix21_long_nom(idx_long),-argo_depth,gvel_NP,0.1:0.1:mv,'LineColor',[0 0 0]+0.5,'LineWidth',1.5) %+ve
c = colorbar;
caxis([-mv mv])
c.Ticks = [-mv:0.1:mv];
c.TickLabels = {'-0.8','','','','-0.4','','','','0','','','','0.4','','','','0.8'}; %mv=0.8
c.TickLength = cticklength;
ylabel(c,'Velocity [m/s]')
box off
colormap(gca,brewermap(16,'*RdBu'))
%labels and axis
text(31.4,-1750,'(b) NP','FontWeight','bold','FontSize',fsize,'Color','w')
XX = [31:1:37];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-1900 -1300 -700 -100];
yticks(YY)
YT = compose('%.0f m',-YY);
yticklabels(YT)
ylabel('Depth [m]')
set(gca,'TickDir','out','FontSize',fsize)


%non-NP cross-section composite
subaxis(3,2,5)
hold on
set(gca,'Color','k') %background colour (and therefore NaNs)
contourf(ix21_long_nom(idx_long),-argo_depth,gvel_normal,-mv:0.1:mv,'LineColor','none') 
contour(ix21_long_nom(idx_long),-argo_depth,gvel_normal,[0 0],'LineColor',[0 0 0],'LineWidth',1.5) %0
contour(ix21_long_nom(idx_long),-argo_depth,gvel_normal,-mv:0.1:-0.1,'--','LineColor',[0 0 0]+0.5,'LineWidth',1.5) %-ve
contour(ix21_long_nom(idx_long),-argo_depth,gvel_normal,0.1:0.1:mv,'LineColor',[0 0 0]+0.5,'LineWidth',1.5) %+ve
c = colorbar;
caxis([-mv mv])
c.Ticks = [-mv:0.1:mv];
c.TickLabels = {'-0.8','','','','-0.4','','','','0','','','','0.4','','','','0.8'}; %mv=0.8
c.TickLength = cticklength;
ylabel(c,'Velocity [m/s]')
box off
colormap(gca,brewermap(16,'*RdBu'))
%labels and axis
text(31.4,-1750,'(c) non-NP','FontWeight','bold','FontSize',fsize,'Color','w')
XX = [31:1:37];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-1900 -1300 -700 -100];
yticks(YY)
YT = compose('%.0f m',-YY);
yticklabels(YT)
ylabel('Depth [m]')
set(gca,'TickDir','out','FontSize',fsize)


%Cross-section composite difference
mv = max(ceil(abs(gvel_diff)*10),[],'all')/10;
subaxis(3,2,4)
hold on
set(gca,'Color','k') %background colour (and therefore NaNs)
contourf(ix21_long_nom(idx_long),-argo_depth,gvel_diff,256,'LineColor','none') 
contour(ix21_long_nom(idx_long),-argo_depth,gvel_diff,[0 0],'LineColor',[0 0 0],'LineWidth',1.5) %0 m/s
contour(ix21_long_nom(idx_long),-argo_depth,gvel_diff,-mv:0.1:-0.1,'--','LineColor',[0 0 0]+0.5,'LineWidth',1.5) %-ve m/s
contour(ix21_long_nom(idx_long),-argo_depth,gvel_diff,0.1:0.1:mv,'LineColor',[0 0 0]+0.5,'LineWidth',1.5) %+ve m/s
colormap(gca,brewermap(256,'*BrBG'))
c = colorbar;
caxis([-mv mv])
c.Ticks = [-mv:0.1:mv];
c.TickLabels = {'-0.4','','-0.2','','0','','0.2','','0.4'}; %mv=0.4
c.TickLength = cticklength;
ylabel(c,'Velocity [m/s]')
box off
%hatching indicates non-significance
[c2,h2] = contour(ix21_long_nom(idx_long),-argo_depth,CI_mask,[1 1]);
set(h2,'linestyle','none','Tag','HatchingRegion');
ax1 = gca;
hp = findobj(ax1,'Tag','HatchingRegion');
hh = hatchfill2(hp,'single','HatchColor','k','HatchDensity',20);
% [X,Y] = meshgrid(ix21_long_nom(idx_long),-argo_depth);
% stipple(X,Y,CI_mask == 1)
%labels and axis
text(31.4,-1750,'(d) Diff.','FontWeight','bold','FontSize',fsize,'Color','w')
XX = [31:1:37];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-1900 -1300 -700 -100];
yticks(YY)
YT = compose('%.0f m',-YY);
yticklabels(YT)
ylabel('Depth [m]')
set(gca,'TickDir','out','FontSize',fsize)


%SLA NP event composite
subaxis(3,2,6)
set(gca,'color',0.6*[1 1 1]) 
hold on
contourf(long_sat_monthly,lat_sat_monthly,mean(sla_monthly0419(:,:,NP_idx),3)',[-0.2:0.025:0.2],'LineColor','none')
contour(long_sat_monthly,lat_sat_monthly,mean(sla_monthly0419(:,:,NP_idx),3)',[0 0],'w','LineWidth',3) %0 m contour 
contour(long_sat_monthly,lat_sat_monthly,mean(sla_monthly0419(:,:,NP_idx),3)',[0.1 0.1],'k--','LineWidth',2) %0.1 m contour
contour(long_sat_monthly,lat_sat_monthly,mean(sla_monthly0419(:,:,NP_idx),3)',[-0.1 -0.1],'k','LineWidth',2) %-0.1 m contour
plot(ix21_long_nom(idx_long),ix21_lat_nom(idx_long),'k','LineWidth',5) %region of IX21 transect shown in composites
plot(ix21_long_nom,ix21_lat_nom,'k','LineWidth',1) %IX21 transect
plot(ix21_mean_core_long,ix21_core_lat,'y.','MarkerSize',25) %mean location of Agulhas core
colormap(gca,brewermap(16,'*PiYG'))
c = colorbar;
caxis([-0.2 0.2])
c.Ticks = [-0.2:0.025:0.2];
c.TickLabels = {'-0.2','','','','-0.1','','','','0','','','','0.1','','','','0.2'};
c.TickLength = cticklength;
c.Label.String = 'SLA [m]';
xlim([25 48])
ylim([-35 -23])
daspect([1 1 1])
box on
%axis ticks and labels:
text(25.9,-24,'(e)','FontWeight','bold','FontSize',fsize,'Color','w')
XX = [25:5:50];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-35:5:-20];
yticks(YY)
YT = compose('%.0f\\circS',-YY);
yticklabels(YT)
set(gca,'TickDir','in','FontSize',fsize)
%grid
for i=1:length(XX)
    xline(XX(i),':')
end
for i=1:length(YY)
    yline(YY(i),':')
end


