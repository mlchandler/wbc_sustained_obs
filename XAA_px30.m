% Mitchell Chandler, SIO
% Last updated: 18/11/2021

%% Load files
load cmap_Zilberman %velocity colourmap
load nom_trans_04-19.mat %nominal transects and coastline

%Bathymetry
topo_long = ncread('etopo6.nc','LON6');
topo_lat = ncread('etopo6.nc','LAT6101_1700');
topo_bath = ncread('etopo6.nc','BATH6');

load LaCasce2020_Ld %first surface mode deformation radius

%% -- Process HR-XBT transects --
%% Read in XBT data
%Read in XBT data
[xbt_p30_long,xbt_p30_lat,xbt_depth,xbt_p30_time,xbt_p30_T] = read_XBT_tem('p30tem9105-1910.nc');
[~,~,~,~,xbt_p30_S] = read_XBT_sal('p30sal9105-1910_argo+argo_corr.nc'); %S

%Some outlier transects that visit different ports at Fiji but as these are
%at the eastern edge I am retaining them for increased data in the WBC

%Restrict to same longitude range as nominal transect
xbt_p30_lat_main = xbt_p30_lat(xbt_p30_long>=p30_long_nom(1) & xbt_p30_long<= p30_long_nom(end),:);
xbt_p30_T_main = xbt_p30_T(xbt_p30_long>=p30_long_nom(1) & xbt_p30_long<= p30_long_nom(end),:,:);
xbt_p30_S_main = xbt_p30_S(xbt_p30_long>=p30_long_nom(1) & xbt_p30_long<= p30_long_nom(end),:,:);
xbt_p30_time_main = xbt_p30_time;

%PX30 is kinked at 168.55

%% Remove transects?
%Remove the 3 transects that start further north
%(24-Nov-2011)
xbt_p30_lat_main(:,find(xbt_p30_time_main == datenum('24-Nov-2011'))) = [];
xbt_p30_T_main(:,:,find(xbt_p30_time_main == datenum('24-Nov-2011'))) = [];
xbt_p30_S_main(:,:,find(xbt_p30_time_main == datenum('24-Nov-2011'))) = [];
xbt_p30_time_main(find(xbt_p30_time_main == datenum('24-Nov-2011'))) = [];

%(03-May-2012)
xbt_p30_lat_main(:,find(xbt_p30_time_main == datenum('03-May-2012'))) = [];
xbt_p30_T_main(:,:,find(xbt_p30_time_main == datenum('03-May-2012'))) = [];
xbt_p30_S_main(:,:,find(xbt_p30_time_main == datenum('03-May-2012'))) = [];
xbt_p30_time_main(find(xbt_p30_time_main == datenum('03-May-2012'))) = [];

%(05-Jun-2013)
xbt_p30_lat_main(:,find(xbt_p30_time_main == datenum('05-Jun-2013'))) = [];
xbt_p30_T_main(:,:,find(xbt_p30_time_main == datenum('05-Jun-2013'))) = [];
xbt_p30_S_main(:,:,find(xbt_p30_time_main == datenum('05-Jun-2013'))) = [];
xbt_p30_time_main(find(xbt_p30_time_main == datenum('05-Jun-2013'))) = [];

%% Plot all XBT T
% for t=1:length(xbt_p30_time_main)
%     figure()
%     contourf(p30_long_nom,-xbt_depth,xbt_p30_T_main(:,:,t)')
%     colorbar
% end

%27-Sep-2005 transect has an unphysical feature about halfway across so
%exclude data eastwards of 165 E
xbt_p30_T_main(p30_long_nom > 165,:,xbt_p30_time_main == datenum('27-Sep-2005')) = NaN;
xbt_p30_S_main(p30_long_nom > 165,:,xbt_p30_time_main == datenum('27-Sep-2005')) = NaN;
xbt_p30_lat_main(p30_long_nom > 165,xbt_p30_time_main == datenum('27-Sep-2005')) = NaN;


%Plot hovmoller of uncorrected XBT surface T to identify where there
%aren't measurements
figure()
contourf(p30_long_nom,xbt_p30_time_main,squeeze(xbt_p30_T_main(:,1,:))','LineColor','none')
datetick('y','keeplimits')
colorbar

%% Plot transects
long_array_use = repmat(p30_long_nom,[1,length(xbt_p30_time_main)]);
long_array_all = repmat(xbt_p30_long,[1,length(xbt_p30_time)]);

%Restrict bathymetry range
topo_bath2 = topo_bath(topo_long>=140 & topo_long<=185,topo_lat>=-35 & topo_lat<=-10)';
topo_long2 = topo_long(topo_long>=140 & topo_long<=185);
topo_lat2 = topo_lat(topo_lat>=-35 & topo_lat<=-10);

%Plot
figure()
hold on
%placeholder plots to get colour key correct
plot(1,1,'r','LineWidth',2)
plot(1,1,'Color','k')
plot(1,1,'Color',rgb('silver'))
%actual plots
contour(topo_long2,topo_lat2,topo_bath2,[-200 -200],'LineColor',rgb('aquamarine')) %200 m isobath
contour(topo_long2,topo_lat2,topo_bath2,[-1000 -1000],'LineColor',[52, 183, 235]/256) %1000 m isobath
contour(topo_long2,topo_lat2,topo_bath2,[-2000 -2000],'LineColor',rgb('blue')) %2000 m isobath
plot(long_array_all,xbt_p30_lat,'Color',rgb('silver'))
plot(long_array_use,xbt_p30_lat_main,'Color','k')
plot(p30_long_nom,p30_lat_nom,'r','LineWidth',2)
plot(gshhs_long,gshhs_lat,'Color',rgb('forest'),'LineWidth',1)
xlim([140 185])
ylim([-30 -10])
daspect([1 1 1])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')
box on
legend('Nominal transect','XBT transects used','Outlier XBT transects','200 m isobath','1000 m isobath','2000 m isobath','Location','SouthOutside')
hold off

%% -- Argo processing steps  --
T = xbt_p30_T_main;
S = xbt_p30_S_main;
xbt_lat = xbt_p30_lat_main;
xbt_lat_nom = p30_lat_nom;
xbt_long_nom = p30_long_nom;
xbt_time = xbt_p30_time_main;
sla_fname = 'sla_daily_2004-2019_px30.nc';
trans_num = 'p30';
argo_fname = 'p30_04-19_delayed.nc';
kink_long = 168.55;
ref_depth = 800;
xx = [140 185];
yy = [-30 -10];

[long_nom,lat_nom,argo_depth,time_monthly,gvel_nom_LNM] = XAA_processing(xbt_lat_nom,xbt_long_nom,xbt_lat,xbt_time,xbt_depth,T,S,sla_fname,trans_num,argo_fname,kink_long,ref_depth,xx,yy);

%% Plot mean geostrophic velocity
gvel_mean_LNM = nanmean(gvel_nom_LNM,3);

%Plot
figure()
mv = ceil(max(abs(gvel_mean_LNM),[],'all')*10)/10;
hold on
set(gca,'Color','k') %background colour (and therefore NaNs)
contourf(long_nom,-argo_depth,gvel_mean_LNM,[-mv:(2*mv)/64:mv],'LineColor','none')
contour(long_nom,-argo_depth,gvel_mean_LNM,[0 0],'LineColor',rgb('silver')) %0 m/s contour
colormap(c_map_Zilberman)
colorbar
caxis([-mv mv])
ylabel(colorbar,'[m/s]')
set(gca,'tickdir','out')
box off

%% -- Reference to Argo sub-surface velocities --
%% Shape of bin for EAC WBC region
figure()
hold on
%nominal transect
plot(long_nom,lat_nom,'r','LineWidth',2,'DisplayName','Nominal transect') 

%trace points on 1000 m isobath
[c_trace] = contour(topo_long2,topo_lat2,topo_bath2,[-1000 -1000],'LineWidth',2,'LineColor',rgb('silver'),'DisplayName','1000-m isobath');
c_trace(:,c_trace(1,:) < 0) = NaN;

%find intersection of nominal transect with 1000-m isobath
[xi,yi] = polyxpoly(long_nom,lat_nom,c_trace(1,:),c_trace(2,:));
plot(xi(1),yi(1),'bo','MarkerFaceColor','b','DisplayName','Transect-Isobath intersection')

%find 1000-m isobath within +/-1.5 lat from intersection point
isobath = c_trace(:,c_trace(2,:) >= yi(1)-1.5 & c_trace(2,:) <= yi(1)+1.5 & c_trace(1,:) < xi(1)+1 & c_trace(1,:) > xi(1)-5);

%find deviation of isobath points from nominal transect intersection
wbc_bath_dx = isobath(1,:) - xi(1);
wbc_bath_dy = isobath(2,:) - yi(1);
plot(xi(1)+wbc_bath_dx,yi(1)+wbc_bath_dy,'g','LineWidth',2,'DisplayName','Isobath segment')

%coastline
plot(gshhs_long,gshhs_lat,'Color',rgb('forest'),'LineWidth',2,'DisplayName','Coastline')

xlim([150 160])
ylim([-30 -20])
daspect([1 1 1])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')
box on
legend('Location','NW')

%% Shape of bin for New Caledonia WBC region
figure()
hold on
%nominal transect
plot(long_nom,lat_nom,'r','LineWidth',2,'DisplayName','Nominal transect') 

%1000-m isobath
plot(c_trace(1,:),c_trace(2,:),'k','DisplayName','1000-m isobath');

%intersection of nominal transect with 1000-m isobath
plot(xi(4),yi(4),'bo','MarkerFaceColor','b','DisplayName','Transect-Isobath intersection')

%trace 1000-m isobath +3 lat along NC east coast
isobath = c_trace(:,c_trace(2,:) >= yi(4) & c_trace(2,:) <= yi(4)+3 & c_trace(1,:) < xi(4) & c_trace(1,:) > xi(4)-5);
%remove trace on other nearby contours
isobath(:,isobath(2,:) > -21.7 & isobath(1,:) > 166.7) = NaN;
isobath(:,isobath(2,:) > -20.8 & isobath(1,:) > 165.6) = NaN;
isobath(:,isobath(1,:) < 164.3) = NaN;
isobath(:,isobath(2,:) < -20.8 & isobath(1,:) < 165.5) = NaN;
isobath(:,isobath(2,:) < -21.7 & isobath(1,:) < 166.74) = NaN;
plot(isobath(1,:),isobath(2,:),'g','LineWidth',2,'DisplayName','Isobath segment')

%linearise segment
p = polyfit(rmmissing(isobath(1,:)),rmmissing(isobath(2,:)),1); %p(1) is gradient p(2) is intercept
lat_hat = polyval(p,rmmissing(isobath(1,:))); %linear fit
plot(rmmissing(isobath(1,:)),lat_hat,'b','LineWidth',2,'DisplayName','Linearised isobath segment')

%find deviation from central point
NC_bath_dx1 = rmmissing(isobath(1,:)) - mean(rmmissing(isobath(1,:)));
NC_bath_dy1 = lat_hat - mean(lat_hat);

%extrapolate linear fit to ensure coverage over +/-1.2 lat range (to reduce the area in the bin)
NC_bath_dy = [-1.2 0 1.2];
NC_bath_dx = interp1(NC_bath_dy1,NC_bath_dx1,NC_bath_dy,'Linear','extrap');

%coastline
plot(gshhs_long,gshhs_lat,'Color',rgb('forest'),'LineWidth',2,'DisplayName','Coastline')

xlim([162 172])
ylim([-25 -17])
daspect([1 1 1])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')
box on
legend('Location','NE')

%% Read in Argo trajectories
[traj_time,traj_long,traj_lat,traj_p,traj_speed,traj_uvel,traj_vvel] = process_argo_traj_v2();

%% Initial quiver plot to identify where the WBC shaped bins are best placed
%(recall PX30 is kinked at 168.55)
kink = find(long_nom>168.5 & long_nom<168.6);
long_nom_box = [xi(1), long_nom(kink), long_nom(end),...
    long_nom(end), long_nom(kink), xi(1),...
    xi(1)];
lat_nom_box = [yi(1)+1.5, lat_nom(kink)+1.5, lat_nom(end)+1.5,...
    lat_nom(end)-1.5, lat_nom(kink)-1.5, yi(1)-1.5,...
    yi(1)+1.5];

%Binning
in = inpolygon(traj_long,traj_lat,long_nom_box,lat_nom_box);
binned_traj_long = traj_long(in);
binned_traj_lat = traj_lat(in);
binned_traj_uvel = traj_uvel(in);
binned_traj_vvel = traj_vvel(in);

figure()
hold on
plot(long_nom,lat_nom,'r')
plot(long_nom_box,lat_nom_box,'k')
quiver(binned_traj_long,binned_traj_lat,binned_traj_uvel,binned_traj_vvel,'b')
plot(gshhs_long,gshhs_lat,'Color',rgb('forest'),'LineWidth',2,'DisplayName','Coastline')
xlim([150 180])
ylim([-30 -15])
daspect([1 1 1])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')
box on

%% Bin trajectories
step = 1/2; %bin width
long_bin = xi(1):step:long_nom(end)+step; %start the binning from the 1000-m isobath (velocity in shallower depths will be referenced to bathymetry LNM)
lat_bin = interp1(long_nom,lat_nom,long_bin,'linear','extrap');

%initialise arrays - each column is a bin
pholder = 1000; %placeholder value
binned_traj_long = NaN(pholder,length(long_bin)-1);
binned_traj_lat = NaN(pholder,length(long_bin)-1);
binned_traj_time = NaN(pholder,length(long_bin)-1);
binned_traj_p = NaN(pholder,length(long_bin)-1);
binned_traj_uvel = NaN(pholder,length(long_bin)-1);
binned_traj_vvel = NaN(pholder,length(long_bin)-1);
mid_long_bin = NaN(1,length(long_bin)-1);
bin_area = NaN(1,length(long_bin)-1);

%initialise plot
figure()
hold on
plot(long_nom,lat_nom,'r','LineWidth',2,'DisplayName','Nominal transect')
contour(topo_long2,topo_lat2,topo_bath2,[-1000 -1000],'LineWidth',2,'LineColor',rgb('silver'),'DisplayName','1000-m isobath')

%binning loop
for i = 1: length(long_bin)-1
    if long_bin(i) < 157 %bins that follow the 1000-m isobath in the EAC region
        xb = [long_bin(i)+wbc_bath_dx long_bin(i+1)+fliplr(wbc_bath_dx) long_bin(i)+wbc_bath_dx(1)];
        yb = [lat_bin(i)+wbc_bath_dy lat_bin(i+1)+fliplr(wbc_bath_dy) lat_bin(i)+wbc_bath_dy(1)];
        
    elseif long_bin(i) < 157+step %match bins between EAC and interior
        xb = [long_bin(i)+wbc_bath_dx long_bin(i+1) long_bin(i+1) long_bin(i)+wbc_bath_dx(1)];
        yb = [lat_bin(i)+wbc_bath_dy lat_bin(i+1)+1.5 lat_bin(i+1)-1.5 lat_bin(i)+wbc_bath_dy(1)];
        
    elseif long_bin(i) >= 166 & long_bin(i) <= 169.5 %bins that follow the 1000-m isobath for New Caledonia
        xb = [long_bin(i)+NC_bath_dx long_bin(i+1)+fliplr(NC_bath_dx) long_bin(i)+NC_bath_dx(1)];
        yb = [lat_bin(i)+NC_bath_dy lat_bin(i+1)+fliplr(NC_bath_dy) lat_bin(i)+NC_bath_dy(1)];
        
    else %interior bins
        xb = [long_bin(i) long_bin(i) long_bin(i+1) long_bin(i+1) long_bin(i)];
        yb = [lat_bin(i)-1.5 lat_bin(i)+1.5 lat_bin(i+1)+1.5 lat_bin(i+1)-1.5 lat_bin(i)-1.5];
    end
    
    %find mid point of each bin
    mid_long_bin(i) = mean([long_bin(i) long_bin(i+1)]);
    
    %binning
    in = inpolygon(traj_long,traj_lat,xb,yb);
    nn = numel(traj_long(in));
    binned_traj_long(1:nn,i) = traj_long(in);
    binned_traj_lat(1:nn,i) = traj_lat(in);
    binned_traj_time(1:nn,i) = traj_time(in);
    binned_traj_p(1:nn,i) = traj_p(in);
    binned_traj_uvel(1:nn,i) = traj_uvel(in);
    binned_traj_vvel(1:nn,i) = traj_vvel(in);
    
    plot(xb,yb,'k') %plot bins
    
    bin_area(i) = polyarea(xb,yb); %computes the area of each bin
end

%find latitude along nominal transect of bin mid-point longs
mid_lat_bin = interp1(long_nom,lat_nom,mid_long_bin);

plot(binned_traj_long,binned_traj_lat,'.') %plot Argo trajectory locations
plot(mid_long_bin,mid_lat_bin,'kx') %plot bin mid-point locations
plot(gshhs_long,gshhs_lat,'Color',rgb('forest'),'LineWidth',2) %plot coastline
xlim([150 180])
ylim([-30 -15])
daspect([1 1 1])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')
box on

%% Binning stats
figure()

%Area of each bin
subplot(2,2,1)
plot(bin_area,'ko')
ylim([1 2])
xlabel('Bin #')
ylabel('Bin area')

%Number of trajectories in each bin
binned_num = NaN(1,size(binned_traj_long,2));
for i=1:size(binned_traj_long,2)
    binned_num(i) = numel(rmmissing(binned_traj_long(:,i)));
end

subplot(2,2,2)
stem(mid_long_bin,binned_num,'k')
ylabel('# trajectories')
xlabel('Longitude [\circE]')
xlim([floor(mid_long_bin(1))-1 ceil(mid_long_bin(end))+1])
set(gca,'YGrid','on','XGrid','off')

%Dates of trajectories in each bin
subplot(2,2,3)
plot(binned_traj_time','kx')
xticks(1:size(binned_traj_time,2))
xlim([0 size(binned_traj_time,2)+1])
datetick('y','keeplimits')
ylabel('Date')
xlabel('Bin #')

%Mean parking depth in each bin
binned_traj_p_mean = nanmean(binned_traj_p);
binned_traj_p_std = nanstd(binned_traj_p);

subplot(2,2,4)
errorbar(mid_long_bin,binned_traj_p_mean,binned_traj_p_std,'k.','MarkerSize',15)
yline(1000,'--')
ylim([950 1050])
xlim([floor(mid_long_bin(1))-1 ceil(mid_long_bin(end))+1])
xlabel('Longitude [\circE]')
ylabel('Depth [dbar]')

%% Velocity vectors with bins overlain
figure()
hold on
quiver(binned_traj_long,binned_traj_lat,binned_traj_uvel,binned_traj_vvel) %plot velocity vectors
plot(gshhs_long,gshhs_lat,'Color',rgb('forest'),'LineWidth',2)
%bins
for i = 1: length(long_bin)-1
    if long_bin(i) < 157 %bins that follow the 1000-m isobath in the EAC region
        xb = [long_bin(i)+wbc_bath_dx long_bin(i+1)+fliplr(wbc_bath_dx) long_bin(i)+wbc_bath_dx(1)];
        yb = [lat_bin(i)+wbc_bath_dy lat_bin(i+1)+fliplr(wbc_bath_dy) lat_bin(i)+wbc_bath_dy(1)];
        
    elseif long_bin(i) < 157+step %match bins between WBC and interior
        xb = [long_bin(i)+wbc_bath_dx long_bin(i+1) long_bin(i+1) long_bin(i)+wbc_bath_dx(1)];
        yb = [lat_bin(i)+wbc_bath_dy lat_bin(i+1)+1.5 lat_bin(i+1)-1.5 lat_bin(i)+wbc_bath_dy(1)];
        
    elseif long_bin(i) >= 166 & long_bin(i) <= 169.5 %bins that follow the 1000-m isobath for New Caledonia
        xb = [long_bin(i)+NC_bath_dx long_bin(i+1)+fliplr(NC_bath_dx) long_bin(i)+NC_bath_dx(1)];
        yb = [lat_bin(i)+NC_bath_dy lat_bin(i+1)+fliplr(NC_bath_dy) lat_bin(i)+NC_bath_dy(1)];
        
    else %interior bins
        xb = [long_bin(i) long_bin(i) long_bin(i+1) long_bin(i+1) long_bin(i)];
        yb = [lat_bin(i)-1.5 lat_bin(i)+1.5 lat_bin(i+1)+1.5 lat_bin(i+1)-1.5 lat_bin(i)-1.5];
    end
    plot(xb,yb,'k') %plot bins
end
contour(topo_long2,topo_lat2,topo_bath2,[-1000 -1000],'LineColor',[1 1 1]*0.7) %1000 m isobath
xlim([150 180])
ylim([-30 -15])
daspect([1 1 1])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')
box on

%% Mean velocity with standard deviation error ellipses
%u
binned_traj_uvel_mean = nanmean(binned_traj_uvel);
%v
binned_traj_vvel_mean = nanmean(binned_traj_vvel);

figure()
hold on
plot(long_nom,lat_nom,'r','LineWidth',1) %plot nominal transect
plot(gshhs_long,gshhs_lat,'Color',rgb('forest'),'LineWidth',2) %plot coastline

sc=10; %scaling factor

%Error ellipses:
%Based heavily on Preisendorfer and Mobley. (1988). Principal Component Analysis in
%Meteorology and Oceanography. Section 21, particularly pages 16-17.
t=-pi:0.01:pi; %for plotting ellipses
for i=1:length(mid_long_bin)
    %remove time-mean for computing (co)variance
    u=rmmissing(binned_traj_uvel(:,i))-nanmean(binned_traj_uvel(:,i));
    v=rmmissing(binned_traj_vvel(:,i))-nanmean(binned_traj_vvel(:,i));
    
    %compute variance and covariance
    sxx = mean(u.*u);
    syy = mean(v.*v);
    sxy = mean(u.*v);
    
    %find principal angle
    thetap = -1/2*atan2(2*sxy,sxx-syy); %radians
    
    %principal variances (along the major and minor axis)
    s11 = 1/2*((sxx+syy) + sqrt((sxx-syy)^2 + 4*sxy^2));
    s22 = 1/2*((sxx+syy) - sqrt((sxx-syy)^2 + 4*sxy^2));
    %take sqrt for std
    std11 = sqrt(s11);
    std22 = sqrt(s22);
    
    %compute ellipse
    xe = std11*cos(t);
    ye = std22*sin(t);
    
    %rotate ellipse by principal angle
    xer = xe*cos(thetap) + ye*sin(thetap);
    yer = -xe*sin(thetap) + ye*cos(thetap);
    
    %centre ellipse around arrow head and apply scaling
    x0=mid_long_bin(i)+sc*binned_traj_uvel_mean(i);
    y0=mid_lat_bin(i)+sc*binned_traj_vvel_mean(i);
    xel = x0+xer*sc;
    yel = y0+yer*sc;
    
    plot(xel,yel,'Color',[0.4 0.4 0.4])
end

%reference ellipse
a=sc*0.1; %x radius
b=sc*0.1; %y radius
x0=175; % x0,y0 centre coordinates
y0=-27;
x=x0+a*cos(t);
y=y0+b*sin(t);
plot(x,y,'Color',[0.4 0.4 0.4]) %plot ellipse
text(176,-27,'u = v = +/- 0.1 m/s','FontSize',12)

%mean velocity vectors
quiver([mid_long_bin 175],[mid_lat_bin -27],sc*[binned_traj_uvel_mean 0],sc*[binned_traj_vvel_mean 0.1],0,'k','LineWidth',2,'ShowArrowHead','on')
text(174.5,-25.5,'0.1 m/s','FontSize',12)

xlim([150 180])
ylim([-29 -18])
daspect([1 1 1])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')
box on

%% Component of Argo velocity perpendicular to nominal transect
%Find the perpendicular component of each individual trajectory and then
%find the mean and SE for each bin 

%initialise mean and SE
vp_mean = NaN*mid_long_bin;
vp_SE = NaN*vp_mean;

%rotate all trajectories and compute mean and SE
for i=1:length(mid_long_bin)
    u = binned_traj_uvel(:,i);
    v = binned_traj_vvel(:,i);
    %find magnitude of vector
    vmag = sqrt(u.^2 + v.^2);
    %find angle of velocity vector with u (horizontal)
    phi = atan2d(v,u);
    
    if mid_long_bin(i) < long_nom(kink) %West of kink
        %find angle of transect with horizontal
        [~,phaseangle] = sw_dist([lat_nom(1) lat_nom(kink)],[long_nom(1) long_nom(kink)]); %E=0, N=90, S=-90
        %subtract 'angle of nominal transect with horizontal' from 'angle of velocity vector from horizontal' to obtain 'angle of velocity vector from nominal transect'
        theta = phi-phaseangle;
        %find component of vector perpendicular to nominal transect
        vp = sind(theta).*vmag;
        %compute mean and SE
        vp_mean(i) = nanmean(vp);
        N = numel(find(isfinite(vp)));
        vp_SE(i) = nanstd(vp)/sqrt(N);
        
    elseif mid_long_bin(i) > long_nom(kink) %East of kink
        %find angle of transect from horizontal
        [~,phaseangle] = sw_dist([lat_nom(kink) lat_nom(end)],[long_nom(kink) long_nom(end)]); %E=0, N=90, S=-90
        %subtract 'angle of nominal transect with horizontal' from 'angle of velocity vector from horizontal' to obtain 'angle of velocity vector from nominal transect'
        theta = phi-phaseangle;
        %find component of vector perpendicular to nominal transect
        vp = sind(theta).*vmag;
        vp_mean(i) = nanmean(vp);
        N = numel(find(isfinite(vp)));
        vp_SE(i) = nanstd(vp)/sqrt(N);
    end
end

%% Plot mean velocity vector and mean perpendicular velocity vector
%calculate velocity components perpendicular to the nominal transect
u = NaN*vp_mean; %perpendicular components then mean
v = NaN*vp_mean;
%first segment
[~,a] = sw_dist([lat_nom(1) lat_nom(kink)],[long_nom(1) long_nom(kink)]);
theta1 = 90-a;
u(mid_long_bin<long_nom(kink)) = -cosd(theta1)*vp_mean(mid_long_bin<long_nom(kink));
v(mid_long_bin<long_nom(kink)) = sind(theta1)*vp_mean(mid_long_bin<long_nom(kink));
%second segment
[~,a] = sw_dist([lat_nom(kink) lat_nom(end)],[long_nom(kink) long_nom(end)]);
theta2 = 90-a;
u(mid_long_bin>long_nom(kink)) = -cosd(theta2)*vp_mean(mid_long_bin>long_nom(kink));
v(mid_long_bin>long_nom(kink)) = sind(theta2)*vp_mean(mid_long_bin>long_nom(kink));

figure() 
hold on
%bathymetry
contour(topo_long2,topo_lat2,topo_bath2,[-1000 -1000],'LineWidth',2,'LineColor',rgb('silver')) %1000 m isobath
plot(gshhs_long,gshhs_lat,'Color',rgb('forest'),'LineWidth',1) %coastline 
plot(long_nom,lat_nom,'Color',[1 1 1]*0.8) %nominal transect
scale = 70;
quiver([mid_long_bin],[mid_lat_bin],scale*[binned_traj_uvel_mean],scale*[binned_traj_vvel_mean],'Color',[1 1 1]*0,'Autoscale','off','ShowArrowHead','off','LineWidth',1) %mean velocity vectors
%velocity perpendicular to nominal transect
pos = v>=0; %northward 
quiver(mid_long_bin(pos),mid_lat_bin(pos),scale*u(pos),scale*v(pos),'Color',[255 160 0]/255,'LineWidth',1,'Autoscale','off','ShowArrowHead','off','LineWidth',2)
neg = v<0; %southward
quiver(mid_long_bin(neg),mid_lat_bin(neg),scale*u(neg),scale*v(neg),'Color',[151 0 255]/255,'LineWidth',1,'Autoscale','off','ShowArrowHead','off','LineWidth',2)
%scale arrow
quiver(175,-27,scale*0.01,scale*0,'Color','k','LineWidth',1,'Autoscale','off','ShowArrowHead','off','LineWidth',1)
text(175,-27.7,'0.01 m/s')
xlim([150 180])
ylim([-30 -18])
daspect([1 1 1])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')
box on
hold off

%% Compare velocity estimates at 1000-m, showing bathymetry along the transect
%compute SE of the Geostrophic velocities at 1000-m (taking into account
%autocorrelation)
gvel_LNM_1000_SE = NaN*long_nom;
for i=1:length(long_nom)
    v = gvel_nom_LNM(argo_depth == 1000,i,:);
    A = xcov(v,'normalized'); %compute autocovariance
    AA = A(round(length(A)/2):end); %consider just the positive half of the record
    max_integral = max(cumtrapz(AA)); %integral timescale approximated by taking the maximum of the autocorrelation integral
    T = 2*max_integral; %integral time scale (multiplied by 2 to take into account +ve and -ve lags)
    eDOF = length(time_monthly)/T;
    if eDOF > length(time_monthly) %maximum DOF is the number of data points
        eDOF = length(time_monthly);
    end
    gvel_LNM_1000_SE(i) = std(v)/sqrt(eDOF);
end

%interpolate bathymetry to nominal transect
bath_nom = interp2(topo_long2,topo_lat2,topo_bath2,long_nom,lat_nom);

%velocities
figure()
subaxis(3,1,[1 2],'SpacingVert',0)
yyaxis left
hold on
errorbar(long_nom,gvel_mean_LNM(argo_depth == 1000,:),gvel_LNM_1000_SE,'LineWidth',1,'DisplayName','Geostrophic velocity (LNM)')
errorbar(mid_long_bin,vp_mean,vp_SE,'-ko','MarkerFaceColor','k','DisplayName','Argo trajectory')
yline(0,'HandleVisibility','off')
xlim([long_nom(1) long_nom(end)])
xticklabels([])
ylabel('Velocity [m/s]')
box on
legend('Location','SE')
set(gca,'YGrid','off','XGrid','on')
yyaxis right
yticks([])

%bathymetry
subaxis(3,1,3)
yyaxis left
yticks([])
yyaxis right
area(long_nom,bath_nom,-6000,'FaceColor','k')
xlim([long_nom(1) long_nom(end)])
yline(-1000,'--')
yline(-2000)
xlabel('Longitude [\circE]')
ylabel('Bathymetry [m]')
set(gca,'YGrid','off','XGrid','on')

%% Re-reference geostrophic velocities to sub-surface velocities
%v(z,x,t) = vg(z,x,t) - <vg(1000,x)> + <vtraj(x)>

%linearly interpolate/extrapolate sub-surface Argo velocities from bin
%mid-points to HR-XBT points 
v_ref = interp1(mid_long_bin,vp_mean,long_nom,'linear','extrap');

v_new = NaN*gvel_nom_LNM; %initialise
for i=1:length(long_nom)
    %don't apply re-referencing where topography is shallower than 1000-m (velocities will be NaN's here already)
    if isnan(mean(gvel_nom_LNM(argo_depth==1000,i,:)))
        v_new(:,i,:) = gvel_nom_LNM(:,i,:);
    else
    %re-reference where topography deeper than 1000-m
    v_new(:,i,:) = gvel_nom_LNM(:,i,:) - mean(gvel_nom_LNM(argo_depth==1000,i,:)) + v_ref(i);
    end
end

gvel_nom = v_new; %save new velocities

gvel_mean = nanmean(gvel_nom,3);

%Plot
figure()
mv = ceil(max(abs(gvel_mean),[],'all')*10)/10;
hold on
set(gca,'Color','k') %background colour (and therefore NaNs)
contourf(long_nom,-argo_depth,gvel_mean,[-mv:(2*mv)/64:mv],'LineColor','none')
contour(long_nom,-argo_depth,gvel_mean,[0 0],'LineColor',rgb('silver')) %0 m/s contour
colormap(c_map_Zilberman)
colorbar
caxis([-mv mv])
ylabel(colorbar,'[m/s]')
set(gca,'tickdir','out')
box off

%% Compare Argo trajectories with re-referenced geostrophic velocity at 1000 m
%compute SE of the LKM geostrophic velocities at 1000-m (taking into account
%autocorrelation)
gvel_1000_SE = NaN*long_nom;
for i=1:length(long_nom)
    v = gvel_nom(argo_depth == 1000,i,:);
    A = xcov(v,'normalized'); %compute autocovariance
    AA = A(round(length(A)/2):end); %consider just the positive half of the record
    max_integral = max(cumtrapz(AA)); %integral timescale approximated by taking the maximum of the autocorrelation integral
    T = 2*max_integral; %integral time scale (multiplied by 2 to take into account +ve and -ve lags)
    eDOF = length(time_monthly)/T;
    if eDOF > length(time_monthly) %maximum DOF is the number of data points
        eDOF = length(time_monthly);
    end
    gvel_1000_SE(i) = std(v)/sqrt(eDOF);
end

figure()
hold on
errorbar(long_nom,gvel_mean_LNM(argo_depth == 1000,:),gvel_LNM_1000_SE,'LineWidth',1,'DisplayName','Geostrophic velocity (LNM)')
errorbar(long_nom,gvel_mean(argo_depth == 1000,:),gvel_1000_SE,'LineWidth',1,'DisplayName','Geostrophic velocity (LKM)')
errorbar(mid_long_bin,vp_mean,vp_SE,'-ko','MarkerFaceColor','k','DisplayName','Argo trajectory')
yline(0,'HandleVisibility','off')
xlim([floor(mid_long_bin(1)) ceil(mid_long_bin(end))])
xlabel('Longitude [\circE]')
ylabel('Velocity [m/s]')
box on
legend('Location','SE')

%% Depth-integrated velocity
%(assuming dbar is equivalent to m)
[depth_int_gvel] = calc_depth_int_vel(long_nom,argo_depth,gvel_nom,time_monthly);
depth_int_gvel_mean = nanmean(depth_int_gvel,2);

%Depth-integrated velocity quiver plot
%calculate velocity components perpendicular to the nominal transect
depth_int_u = NaN*depth_int_gvel_mean;
depth_int_v = NaN*depth_int_gvel_mean;
%first segment
[~,a] = sw_dist([lat_nom(1) lat_nom(find(long_nom<kink_long-0.05,1,'last'))],[long_nom(1) long_nom(find(long_nom<kink_long-0.05,1,'last'))]);
theta1 = 90-a;
depth_int_u(find(long_nom<kink_long-0.05)) = -cosd(theta1)*depth_int_gvel_mean(find(long_nom<kink_long-0.05));
depth_int_v(find(long_nom<kink_long-0.05)) = sind(theta1)*depth_int_gvel_mean(find(long_nom<kink_long-0.05));
%second segment
[~,a] = sw_dist([lat_nom(find(long_nom>kink_long+0.05,1,'first')) lat_nom(end)],[long_nom(find(long_nom>kink_long+0.05,1,'first')) long_nom(end)]);
theta2 = 90-a;
depth_int_u(find(long_nom>kink_long+0.05)) = -cosd(theta2)*depth_int_gvel_mean(find(long_nom>kink_long+0.05));
depth_int_v(find(long_nom>kink_long+0.05)) = sind(theta2)*depth_int_gvel_mean(find(long_nom>kink_long+0.05));
%kink point
depth_int_u(find(long_nom>kink_long-0.05 & long_nom<kink_long+0.05)) = mean([-cosd(theta1)*depth_int_gvel_mean(find(long_nom>kink_long-0.05 & long_nom<kink_long+0.05)), -cosd(theta2)*depth_int_gvel_mean(find(long_nom>kink_long-0.05 & long_nom<kink_long+0.05))]);
depth_int_v(find(long_nom>kink_long-0.05 & long_nom<kink_long+0.05)) = mean([sind(theta1)*depth_int_gvel_mean(find(long_nom>kink_long-0.05 & long_nom<kink_long+0.05)), sind(theta2)*depth_int_gvel_mean(find(long_nom>kink_long-0.05 & long_nom<kink_long+0.05))]);


figure() 
hold on
%bathymetry
contour(topo_long2,topo_lat2,topo_bath2,[-200 -200],'LineColor',rgb('aquamarine')) %200 m isobath
contour(topo_long2,topo_lat2,topo_bath2,[-1000 -1000],'LineColor',[52, 183, 235]/256) %1000 m isobath
contour(topo_long2,topo_lat2,topo_bath2,[-2000 -2000],'LineColor',rgb('blue')) %2000 m isobath
plot(gshhs_long,gshhs_lat,'Color',rgb('forest'),'LineWidth',1,'HandleVisibility','off')
%depth-integrated velocities
scale = 0.02;
% quiver(long_nom,lat_nom,scale*depth_int_u,scale*depth_int_v,'LineWidth',1.5,'Autoscale','off') %to confirm the quiver arrows are all scaled the same
%N-ward vectors
pos = depth_int_v>=0; 
quiver(long_nom(pos),lat_nom(pos),scale*depth_int_u(pos),scale*depth_int_v(pos),'Color',[255 160 0]/255,'LineWidth',1,'Autoscale','off','ShowArrowHead','off')
%S-ward vectors
neg = depth_int_v<0; 
quiver(long_nom(neg),lat_nom(neg),scale*depth_int_u(neg),scale*depth_int_v(neg),'Color',[151 0 255]/255,'LineWidth',1,'Autoscale','off','ShowArrowHead','off')
%Scale arrow
quiver(147,-29,scale*0',scale*100','Color','k','LineWidth',1,'Autoscale','off','ShowArrowHead','on')
text(147.5,-28,'100 m^2/s')
xlim([146 180])
ylim([-33 -16])
daspect([1 1 1])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')
box on
legend('200 m isobath','1000 m isobath','2000 m isobath','Location','SE')
hold off


%% -- Explore WBC indices --
figure() 
mv = round(max(abs(depth_int_gvel),[],'all'),-2); %round to the nearest 100
hold on
contourf(long_nom,time_monthly,depth_int_gvel',[-mv:(2*mv)/64:mv],'Linecolor','none') 
contour(long_nom,time_monthly,depth_int_gvel',[0 0],'LineColor',[0.8 0.8 0.8],'LineWidth',0.5) 
plot(long_nom(end),xbt_time,'k<','MarkerSize',2,'MarkerFaceColor','k') %timestamps for XBT transects
colormap(c_map_Zilberman)
colorbar
caxis([-mv mv])
set(gca,'Color','k') 
yticks(datenum(2005,1:24:204,1))
datetick('y','keeplimits','keepticks')
xlabel('Longitude [\circE]')
ylabel('Year')
ylabel(colorbar,'Depth-integrated geostrophic velocity [m^2/s]')

%% Estimate appropriate distance to smooth over
figure()
contourf(WOA18_Ld.xt,WOA18_Ld.yt,WOA18_Ld.Ld'/1E3,20,'LineColor','none')
colorbar
ylabel(colorbar,'Deformation radius [km]')
hold on
plot(long_nom,lat_nom,'r','LineWidth',2)
plot(gshhs_long,gshhs_lat,'Color',rgb('forest'),'LineWidth',1)
axis equal
xlim([140 185])
ylim([-30 -10])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')


%deformation radius along nominal transect
LD = interp2(WOA18_Ld.xt,WOA18_Ld.yt,WOA18_Ld.Ld',long_nom,lat_nom);
figure()
plot(long_nom,LD/1E3)
xlabel('Longitude [\circE]')
ylabel('Deformation radius [km]')

%% Smooth in longitude
w_size = 5; %window size

W = (1/sum(triang(w_size)))*triang(w_size);
smth_depth_int_gvel=NaN*depth_int_gvel;
for t=1:length(time_monthly)
smth_depth_int_gvel(:,t) = conv_filt(depth_int_gvel(:,t),W,w_size);
end

mv = round(max(abs(depth_int_gvel),[],'all'),-2); %round to the nearest 100

figure() 
hold on
contourf(long_nom,time_monthly,smth_depth_int_gvel',[-mv:(2*mv)/64:mv],'Linecolor','none') 
contour(long_nom,time_monthly,smth_depth_int_gvel',[0 0],'LineColor',[0.8 0.8 0.8],'LineWidth',0.5) 
plot(long_nom(end),xbt_time,'k<','MarkerSize',4,'MarkerFaceColor','w') %timestamps for XBT transects
colormap(c_map_Zilberman)
colorbar
caxis([-mv mv])
set(gca,'Color','k') 
yticks(datenum(2005,1:24:204,1))
datetick('y','keeplimits','keepticks')
xlabel('Longitude [\circE]')
ylabel('Year')
ylabel(colorbar,'Depth-integrated geostrophic velocity [m^2/s]')

%% Identify core and inshore and offshore edges of WBC
%Find location of jet core
[core_speed,core_idx] = min(smth_depth_int_gvel,[],1);
core_long = long_nom(core_idx);
%mean latitude corresponding to the mean jet core longitude
mean_core_lat = round(interp1(long_nom,lat_nom,mean(core_long)),1); 

%Find offshore edge of the WBC
offshore_long = NaN(length(time_monthly),1);
offshore_idx = offshore_long;
for t=1:length(time_monthly)
    count=1;
    zerox_idx = find(smth_depth_int_gvel(1:end-1,t).* smth_depth_int_gvel(2:end,t) <= 0); %find indices before sign change
    %offshore point has to be east of WBC core
    offshore_long(t) = long_nom(zerox_idx(count));
    while offshore_long(t) < core_long(t)
        count=count+1;
        offshore_long(t) = long_nom(zerox_idx(count));
    end
    offshore_idx(t) = zerox_idx(count);
end
%plot the offshore edge longitude half a space further offshore as the
%indices are for before where the sign changes
offshore_long_corrected = offshore_long+mode(diff(long_nom))/2; 

%Find inshore edge of the WBC
%inshore edge is at the western edge of the transect
inshore_long = NaN(length(time_monthly),1);
inshore_idx = inshore_long;
for t=1:length(time_monthly)
        ii = find(~isnan(smth_depth_int_gvel(:,t)),1,'first'); %first non-NaN depth-int gvel value
        inshore_long(t) = long_nom(ii);
     inshore_idx(t) = ii;
end

figure() 
hold on
contourf(long_nom,time_monthly,smth_depth_int_gvel',[-mv:(2*mv)/64:mv],'Linecolor','none') 
contour(long_nom,time_monthly,smth_depth_int_gvel',[0 0],'LineColor',[0.8 0.8 0.8],'LineWidth',0.5) 
plot(long_nom(end),xbt_time,'k<','MarkerSize',4,'MarkerFaceColor','w') %timestamps for XBT transects
plot(core_long,time_monthly,'r--') %location of jet core
plot(offshore_long_corrected,time_monthly,'r-','LineWidth',1.5) %offshore edge
plot(inshore_long,time_monthly,'r-','LineWidth',1) %inshore edge
colormap(c_map_Zilberman)
colorbar
caxis([-mv mv])
set(gca,'Color','k') 
yticks(datenum(2005,1:24:204,1))
datetick('y','keeplimits','keepticks')
xlabel('Longitude [\circE]')
ylabel('Year')
ylabel(colorbar,'Depth-integrated geostrophic velocity [m^2/s]')

%% Save variables for plotting a combined figure
% px30_long_nom = long_nom;
% px30_divel = smth_depth_int_gvel;
% px30_xbt_times = xbt_time;
% px30_core_long = core_long;
% px30_offshore_long = offshore_long_corrected;
% save('px30_fig2.mat',...
%     'time_monthly','px30_long_nom','px30_divel','px30_xbt_times','px30_core_long','px30_offshore_long');

%% Cumulative transport
%using smoothed depth-integrated velocity
[cumulated_transport,transport_long,~] = calc_cumulative_transport(long_nom,lat_nom,smth_depth_int_gvel,time_monthly);
cum_transport_mean = nanmean(cumulated_transport,2);

figure()
plot(transport_long,cumulated_transport)
hold on
plot(transport_long,cum_transport_mean,'k','LineWidth',2)
yline(0,'--')
xlim([transport_long(1) transport_long(end)])
xlabel('Longitude [\circE]')
ylabel('Eastward-accumulated transport [Sv]')

%% Derived indices - transport and deviations in offshore edge
%WBC transport between inshore and offshore edge
wbc_transport = NaN*time_monthly;
for t=1:length(time_monthly)
    wbc_transport(t) = cumulated_transport(offshore_idx(t)+1,t) - cumulated_transport(inshore_idx(t),t);
end

%Compute deviations in distance from the location of the offshore edge
%mean location
offshore_long_c_mean = mean(offshore_long_corrected);
offshore_lat_c_mean = interp1(long_nom,lat_nom,offshore_long_c_mean);

offshore_lat_corrected = interp1(long_nom,lat_nom,offshore_long_corrected);

offshore_deviations = NaN*offshore_long;
for t=1:length(time_monthly)
    offshore_deviations(t) = gsw_distance([offshore_long_c_mean, offshore_long_corrected(t)],[offshore_lat_c_mean, offshore_lat_corrected(t)])/1E3; %km
    %make distance negative for locations inshore (west) of the mean
    if offshore_long_corrected(t) < offshore_long_c_mean
        offshore_deviations(t) = -offshore_deviations(t);
    end
end

%% Plot time series
%Linear trends - monthly time step (dt=1), 95% CI (alpha=0.05).
%trends are not significant if the gradient includes 0 within its CI
%-wbc transport
[p_transport,trend_transport,CI_transport] = linear_trend(time_monthly,wbc_transport,1,0.05);

%-core speed
[p_core,trend_core,CI_core] = linear_trend(time_monthly,core_speed,1,0.05); %have to transpose this - but work out why so I can include in function

%-offshore edge deviations
[p_offshore,trend_offshore,CI_offshore] = linear_trend(time_monthly,offshore_deviations,1,0.05);

%wbc transport
figure()
hold on
plot(time_monthly,wbc_transport,'k','LineWidth',1)
%only plot trend if significant
if abs(p_transport(1)) - abs(CI_transport) > 0 
    plot(time_monthly,trend_transport,'r','LineWidth',1)
end
%print trend and CI values on plot
text(datenum('01-Jun-2004'),-38,{[num2str(round(p_transport(1)*365.25,2)),' +/- ',num2str(round(CI_transport*365.25,2)),' Sv/yr']})
datetick('x')
ylabel('WBC transport [Sv]')
box on

%core speed
figure()
hold on
plot(time_monthly,core_speed,'k','LineWidth',1)
%only plot trend if significant
if abs(p_core(1)) - abs(CI_core) > 0 
    plot(time_monthly,trend_core,'r','LineWidth',1)
end
%print trend and CI values on plot
text(datenum('01-Jun-2004'),-325,{[num2str(round(p_core(1)*365.25,2)),' +/- ',num2str(round(CI_core*365.25,2)),' m^2/s/yr']})
datetick('x')
ylabel('Core depth-integrated velocity [m^2/s]')
box on

%deviation in offshore edge (distance)
figure()
hold on
plot(time_monthly,offshore_deviations,'k','LineWidth',1)
%only plot trend if significant
if abs(p_offshore(1)) - abs(CI_offshore) > 0 
    plot(time_monthly,trend_offshore,'r','LineWidth',1)
end
%print trend and CI values on plot
text(datenum('01-Jun-2004'),275,{[num2str(round(p_offshore(1)*365.25,2)),' +/- ',num2str(round(CI_offshore*365.25,2)),' km/yr']})
yline(0,'--')
datetick('x')
ylabel('Deviations in offshore edge [km]')
box on

%% Annual cycle
%3-month filter of time series
w_size = 3; %window size

%window
% W = (1/sum(triang(w_size)))*triang(w_size); %for window size 3 a hanning filter is almost the same as a triangle filter
W = (1/sum(boxcar(w_size)))*boxcar(w_size); %boxcar is the same as using: movmean(wbc_transport,w_size,'EndPoints','fill');

%filter
%-wbc transport
[wbc_transport_filt3] = conv_filt(wbc_transport,W,w_size);
%-core speed
[core_speed_filt3] = conv_filt(core_speed,W,w_size);
%-offshore edge deviations
[offshore_deviations_filt3] = conv_filt(offshore_deviations,W,w_size);

%Average for each month
[~,mm,~] = ymd(datetime(time_monthly,'ConvertFrom','datenum'));

%initialise 
wbc_transport_month = NaN(12,1);
core_speed_month = NaN(12,1);
offshore_deviations_month = NaN(12,1);

wbc_transport_month_SE = wbc_transport_month;
core_speed_month_SE = core_speed_month;
offshore_deviations_month_SE = offshore_deviations_month;

for i=1:12
%     %no smoothing
%     wbc_transport_month(i) = mean(wbc_transport(mm == i));
%     core_speed_month(i) = mean(core_speed(mm == i));
%     offshore_deviations_month(i) = mean(offshore_deviations(mm == i));
%     
%     %SE = std/sqrt(N)
%     wbc_transport_month_SE(i) = std(wbc_transport(mm == i))/sqrt(length(find(mm == i)));
%     core_speed_month_SE(i) = std(core_speed(mm == i))/sqrt(length(find(mm == i)));
%     offshore_deviations_month_SE(i) = std(offshore_deviations(mm == i))/sqrt(length(find(mm == i)));
    
    
    %with smoothing
    wbc_transport_month(i) = nanmean(wbc_transport_filt3(mm == i));
    core_speed_month(i) = nanmean(core_speed_filt3(mm == i));
    offshore_deviations_month(i) = nanmean(offshore_deviations_filt3(mm == i));
   
    %SE = std/sqrt(N)
    wbc_transport_month_SE(i) = nanstd(wbc_transport_filt3(mm == i))/sqrt(length(find(~isnan(wbc_transport_filt3(mm == i)))));
    core_speed_month_SE(i) = nanstd(core_speed_filt3(mm == i))/sqrt(length(find(~isnan(core_speed_filt3(mm == i)))));
    offshore_deviations_month_SE(i) = nanstd(offshore_deviations_filt3(mm == i))/sqrt(length(find(~isnan(offshore_deviations_filt3(mm == i)))));
end

%Plot
%-wbc transport
figure()
hold on
shaded_error(1:12,wbc_transport_month,wbc_transport_month_SE,'k',0.1,2)
yline(mean(wbc_transport_month),'--')
ylabel('WBC transport [Sv]')
xlim([1 12])
xticks([1:12])
xticklabels(['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'])
box on

%-core speed
figure()
hold on
shaded_error(1:12,core_speed_month,core_speed_month_SE,'k',0.1,2)
yline(mean(core_speed_month),'--')
ylabel('Core depth-integrated velocity [m^2/s]')
xlim([1 12])
xticks([1:12])
xticklabels(['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'])
box on

%-offshore edge deviations
figure()
hold on
shaded_error(1:12,offshore_deviations_month,offshore_deviations_month_SE,'k',0.1,2)
yline(mean(offshore_deviations_month),'--')
ylabel('Deviations in offshore edge [km]')
xlim([1 12])
xticks([1:12])
xticklabels(['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'])
box on

%% Assess significance of annual cycle
%Assess significance of annual cycles by correlating with the 3-month
%smoothed and detrended time series

%Initialise table
ann_cyc_sig = array2table(NaN(3,4),... 
    'VariableNames',{'r','p','r2','R2'},...
    'RowNames',{'Transport','Core speed','Offshore'}); 

%-wbc transport
ann_cyc = detrend(wbc_transport_month(mm));
filtered_ts = detrend(wbc_transport_filt3,'omitnan');
[r,p] = corr_pval(rmmissing(filtered_ts),ann_cyc(2:end-1));
r2 = r^2;
R2 = 1 - nansum((filtered_ts - ann_cyc).^2)/nansum((filtered_ts - nanmean(filtered_ts)).^2); %check that R2 is the same
ann_cyc_sig{'Transport','r'} = r;
ann_cyc_sig{'Transport','p'} = p;
ann_cyc_sig{'Transport','r2'} = r2;
ann_cyc_sig{'Transport','R2'} = R2;

%-core speed
ann_cyc = detrend(core_speed_month(mm));
filtered_ts = detrend(core_speed_filt3','omitnan');
[r,p] = corr_pval(rmmissing(filtered_ts),ann_cyc(2:end-1));
r2 = r^2;
R2 = 1 - nansum((filtered_ts - ann_cyc).^2)/nansum((filtered_ts - nanmean(filtered_ts)).^2); %check that R2 is the same
ann_cyc_sig{'Core speed','r'} = r;
ann_cyc_sig{'Core speed','p'} = p;
ann_cyc_sig{'Core speed','r2'} = r2;
ann_cyc_sig{'Core speed','R2'} = R2;

%-offshore edge deviations
ann_cyc = detrend(offshore_deviations_month(mm));
filtered_ts = detrend(offshore_deviations_filt3,'omitnan');
[r,p] = corr_pval(rmmissing(filtered_ts),ann_cyc(2:end-1));
r2 = r^2;
R2 = 1 - nansum((filtered_ts - ann_cyc).^2)/nansum((filtered_ts - nanmean(filtered_ts)).^2); %check that R2 is the same
ann_cyc_sig{'Offshore','r'} = r;
ann_cyc_sig{'Offshore','p'} = p;
ann_cyc_sig{'Offshore','r2'} = r2;
ann_cyc_sig{'Offshore','R2'} = R2;

%% Time-mean transport
[mean_transport,transport_uncertainty,transport_stdev,transport_SE] =...
    mean_error_transport(wbc_transport,vp_SE,mid_long_bin,long_nom,lat_nom,inshore_long(1),offshore_long_c_mean)

%% -- Save variables --
%% Save transect 
% px30_long_nom = long_nom;
% px30_lat_nom = lat_nom;
% px30_gvel_LKM = gvel_nom;
% 
% save('px30_velocity.mat',...
%     'time_monthly','argo_depth',...
%     'px30_long_nom','px30_lat_nom','px30_gvel_LKM');

%% Save WBC metrics
% px30_wbc_transport_raw = wbc_transport;
% px30_core_speed_raw = core_speed;
% px30_offshore_dev_raw = offshore_deviations;
% 
% px30_wbc_transport_month = wbc_transport_month;
% px30_wbc_transport_month_SE = wbc_transport_month_SE;
% px30_core_speed_month = core_speed_month;
% px30_core_speed_month_SE = core_speed_month_SE;
% px30_offshore_deviations_month = offshore_deviations_month;
% px30_offshore_deviations_month_SE = offshore_deviations_month_SE;
% 
% px30_core_lat = mean_core_lat;
% 
% save('px30_variability.mat',...
%     'time_monthly','px30_wbc_transport_raw','px30_core_speed_raw','px30_offshore_dev_raw',...
%     'px30_wbc_transport_month','px30_wbc_transport_month_SE',...
%     'px30_core_speed_month','px30_core_speed_month_SE',...
%     'px30_offshore_deviations_month','px30_offshore_deviations_month_SE',...
%     'px30_core_lat');

%% Create netCDF for velocity data
% delete px30_velocity.nc
% 
% %Write netcdf variables:
% %time
% nccreate('px30_velocity.nc','time','Dimensions',{'time',192})
% ncwrite('px30_velocity.nc','time',time_monthly-datenum('01-Jan-1970'))
% ncwriteatt('px30_velocity.nc','time','standard_name','time')
% ncwriteatt('px30_velocity.nc','time','long_name','time')
% ncwriteatt('px30_velocity.nc','time','units','days since 1970-01-01 00:00:00')
% ncwriteatt('px30_velocity.nc','time','comments','given as the 15th of each month')
% 
% %nominal transect longitude
% nccreate('px30_velocity.nc','longitude','Dimensions',{'longitude',236})
% ncwrite('px30_velocity.nc','longitude',long_nom)
% ncwriteatt('px30_velocity.nc','longitude','standard_name','longitude')
% ncwriteatt('px30_velocity.nc','longitude','long_name','nominal transect longitude')
% ncwriteatt('px30_velocity.nc','longitude','units','deg_E')
% 
% %nominal transect latitude
% nccreate('px30_velocity.nc','latitude','Dimensions',{'longitude',236})
% ncwrite('px30_velocity.nc','latitude',lat_nom)
% ncwriteatt('px30_velocity.nc','latitude','standard_name','latitude')
% ncwriteatt('px30_velocity.nc','latitude','long_name','nominal transect latitude')
% ncwriteatt('px30_velocity.nc','latitude','units','deg_N')
% 
% %depth
% nccreate('px30_velocity.nc','depth','Dimensions',{'depth',76})
% ncwrite('px30_velocity.nc','depth',argo_depth)
% ncwriteatt('px30_velocity.nc','depth','standard_name','depth')
% ncwriteatt('px30_velocity.nc','depth','long_name','depth')
% ncwriteatt('px30_velocity.nc','depth','units','m')
% ncwriteatt('px30_velocity.nc','depth','positive','down')
% 
% %velocity cross-sections (LKM)
% nccreate('px30_velocity.nc','vel','Dimensions',{'depth',76,'longitude',236,'time',192},'FillValue',NaN)
% ncwrite('px30_velocity.nc','vel',gvel_nom)
% ncwriteatt('px30_velocity.nc','vel','long_name','monthly cross-transect absolute geostrophic velocity')
% ncwriteatt('px30_velocity.nc','vel','units','m/s')
% ncwriteatt('px30_velocity.nc','vel','positive','northward')
% ncwriteatt('px30_velocity.nc','vel','comments','cross-transect geostrophic velocity referenced to Argo trajectory velocity at 1000-m')
% 
% %velocity cross-sections (LNM)
% nccreate('px30_velocity.nc','gvel_LNM','Dimensions',{'depth',76,'longitude',236,'time',192},'FillValue',NaN)
% ncwrite('px30_velocity.nc','gvel_LNM',gvel_nom_LNM)
% ncwriteatt('px30_velocity.nc','gvel_LNM','long_name','monthly cross-transect geostrophic velocity')
% ncwriteatt('px30_velocity.nc','gvel_LNM','units','m/s')
% ncwriteatt('px30_velocity.nc','gvel_LNM','positive','northward')
% ncwriteatt('px30_velocity.nc','gvel_LNM','comments','cross-transect geostrophic velocity referenced to a level-of-no-motion at 1975-m, or the bathymetry if shallower')
% 
% %longitude for Argo trajectory bins
% nccreate('px30_velocity.nc','long_for_vel_err','Dimensions',{'long_for_vel_err',47})
% ncwrite('px30_velocity.nc','long_for_vel_err',mid_long_bin')
% ncwriteatt('px30_velocity.nc','long_for_vel_err','standard_name','longitude')
% ncwriteatt('px30_velocity.nc','long_for_vel_err','long_name','longitude along nominal transect for velocity uncertainty term')
% ncwriteatt('px30_velocity.nc','long_for_vel_err','units','deg_E')
% ncwriteatt('px30_velocity.nc','long_for_vel_err','comments','longitude along the nominal transect at the mid point of each 1/2 degree longitude x 3 degree latitude bin used to compute Argo trajectory velocity at 1000-m')
% 
% %latitude for Argo trajectory bins
% nccreate('px30_velocity.nc','lat_for_vel_err','Dimensions',{'long_for_vel_err',47})
% ncwrite('px30_velocity.nc','lat_for_vel_err',mid_lat_bin')
% ncwriteatt('px30_velocity.nc','lat_for_vel_err','standard_name','latitude')
% ncwriteatt('px30_velocity.nc','lat_for_vel_err','long_name','latitude along nominal transect for velocity uncertainty term')
% ncwriteatt('px30_velocity.nc','lat_for_vel_err','units','deg_N')
% ncwriteatt('px30_velocity.nc','lat_for_vel_err','comments','latitude along the nominal transect at the mid point of each 1/2 degree longitude x 3 degree latitude bin used to compute Argo trajectory velocity at 1000-m')
% 
% %argo velocity standard error for each bin
% nccreate('px30_velocity.nc','vel_err','Dimensions',{'long_for_vel_err',47})
% ncwrite('px30_velocity.nc','vel_err',vp_SE')
% ncwriteatt('px30_velocity.nc','vel_err','long_name','velocity uncertainty term')
% ncwriteatt('px30_velocity.nc','vel_err','units','m/s')
% ncwriteatt('px30_velocity.nc','vel_err','positive','northward')
% ncwriteatt('px30_velocity.nc','vel_err','comments','standard error of Argo trajectory velocity at 1000-m perpendicular to the nominal transect in each 1/2 degree longitude x 3 degree latitude bin')
% 
% 
% %Global attributes:
% ncwriteatt('px30_velocity.nc','/','title','px30_velocity.nc');
% ncwriteatt('px30_velocity.nc','/','summary','cross-sectional time series of absolute geostrophic velocity across HR-XBT line PX30 between Brisbane, Australia and Suva, Fiji from 0-m to 1975-m'); 
% ncwriteatt('px30_velocity.nc','/','time_coverage_start','15-Jan-2004');
% ncwriteatt('px30_velocity.nc','/','time_coverage_end','15-Dec-2019');
% ncwriteatt('px30_velocity.nc','/','latitude_min',min(lat_nom));
% ncwriteatt('px30_velocity.nc','/','latitude_max',max(lat_nom));
% ncwriteatt('px30_velocity.nc','/','longitude_min',min(long_nom));
% ncwriteatt('px30_velocity.nc','/','longitude_max',max(long_nom));
% ncwriteatt('px30_velocity.nc','/','longitude_res','0.1 deg');
% ncwriteatt('px30_velocity.nc','/','depth_min',min(argo_depth));
% ncwriteatt('px30_velocity.nc','/','depth_max',max(argo_depth));
% ncwriteatt('px30_velocity.nc','/','reference','[paper to be added upon acceptance]'); %<<add paper here
% ncwriteatt('px30_velocity.nc','/','date_created',datestr(now));
% ncwriteatt('px30_velocity.nc','/','creator_name','Mitchell Chandler');
% ncwriteatt('px30_velocity.nc','/','creator_email','mlchandl@ucsd.edu');
% ncwriteatt('px30_velocity.nc','/','institution','Scripps Institution of Oceanography');
% ncwriteatt('px30_velocity.nc','/','version','1');
% ncwriteatt('px30_velocity.nc','/','history','version 1 created 18-Nov-2021');
% 
% 
% %Display netcdf:
% ncdisp('px30_velocity.nc')

%% -- End -- 
