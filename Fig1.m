% Mitchell Chandler, SIO
% Last updated: 01/07/2021

load ix21_velocity
load px30_velocity
load px40_velocity

%High-resolution Argo climatology SH
load argo_mean2018_DH_hr.mat %can jump straight to plotting after loading this

%% Compute Argo DH
%Read in high-res Argo climatology 
argo_mean2018_lat = ncread('RG_ArgoClim_33pfit_2019_mean.nc','LATITUDE');
argo_mean2018_long = ncread('RG_ArgoClim_33pfit_2019_mean.nc','LONGITUDE');
argo_mean2018_p = ncread('RG_ArgoClim_33pfit_2019_mean.nc','PRESSURE');
argo_mean2018_T = ncread('RG_ArgoClim_33pfit_2019_mean.nc','ARGO_TEMPERATURE_MEAN');
argo_mean2018_S = ncread('RG_ArgoClim_33pfit_2019_mean.nc','ARGO_SALINITY_MEAN');
%high-resolution SH has also been saved so that it can just be read in and
%plotted because it takes such a long time to run

%Restrict range to Indo-Pacifc
lat_idx = find(argo_mean2018_lat <= 50 & argo_mean2018_lat >= -50);
long_idx = find(argo_mean2018_long <= 300 & argo_mean2018_long >= 20);

argo_T = argo_mean2018_T(long_idx,lat_idx,:);
argo_S = argo_mean2018_S(long_idx,lat_idx,:);
argo_LA = argo_mean2018_lat(lat_idx);
argo_LO = argo_mean2018_long(long_idx);

%Compute DH
%reshape matrices into vectors
dim = size(argo_T);
T_vec = reshape(argo_T,[dim(1)*dim(2),dim(3)]); 
S_vec = reshape(argo_S,[dim(1)*dim(2),dim(3)]); 
Long_vec = repmat(argo_LO,dim(2),1);
argo_LA_rep = repmat(argo_LA,1,dim(1))';
Lat_vec = reshape(argo_LA_rep,[dim(1)*dim(2),1]);

%set NaNs throughout water column for values which are NaN at the bottom
T_vec(isnan(T_vec(:,end)),:) = NaN;
S_vec(isnan(S_vec(:,end)),:) = NaN;

%convert practical salinity to absolute salinity
SA = gsw_SA_from_SP(S_vec,argo_mean2018_p,Long_vec,Lat_vec);
%convert in-situ temperature to conservative temperature
CT = gsw_CT_from_t(SA,T_vec,argo_mean2018_p);

%calculate dynamic height relative to bottom argo p level (1975 dbar)
argo_DH_vec = gsw_geo_strf_dyn_height(SA',CT',argo_mean2018_p,argo_mean2018_p(end));

%reshape DH back into matrix
argo_DH_vec = argo_DH_vec';
argo_DH = reshape(argo_DH_vec,[dim(1),dim(2),dim(3)]);

%% Plot
figure('color','w');
% clf
hold on
%map projection
m_proj('Miller','lon',[20 300],'lat',[-45 45]);
%DH
m_contourf(argo_LO,argo_LA,argo_DH(:,:,1)',20,'LineColor','none')
colormap(flipud(crameri('roma')))
c = colorbar('eastoutside');
c.FontSize = 12;
ylabel(c,'Dynamic Height [m^2 s^{-2}]','FontSize',13);
%coastline
m_coast('patch',[.7 .7 .7],'edgecolor','none');
%nominal transects
m_plot(px30_long_nom,px30_lat_nom,'Color','k','LineWidth',5)
m_plot(px40_long_nom,px40_lat_nom,'Color','k','LineWidth',5)
m_plot(ix21_long_nom,ix21_lat_nom,'Color','k','LineWidth',5)
%boxes identifying regions of other panels
m_plot([30 60 60 30 30],[-21 -21 -36 -36 -21],'k') %IX21 
m_plot([150 180 180 150 150],[-18 -18 -33 -33 -18],'k') %PX30
m_plot([135 165 165 135 135],[42 42 27 27 42],'k') %PX40
%grid
m_grid('linestyle','none','box','fancy','xtick',14,'ytick',7,'tickdir','in','ticklength',2.5E-3,'FontSize',12);

%% 
% save('argo_mean2018_DH_hr.mat','argo_DH','argo_LO','argo_LA','argo_mean2018_p'); %save high-res dynamic height because it takes so long to run

