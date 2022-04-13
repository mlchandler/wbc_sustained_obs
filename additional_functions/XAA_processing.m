% Mitchell Chandler, SIO
% Last updated: 12/04/2022

function [long_nom,lat_nom,argo_p_1975,time_monthly,gvel_nom] = XAA_processing(xbt_lat_nom,xbt_long_nom,xbt_coord,xbt_time,xbt_depth,T,S,sla_fname,trans_num,argo_fname,kink_long,ref_depth,xx,yy)
%% Make sure variables are double format
xbt_lat_nom = double(xbt_lat_nom);
xbt_long_nom = double(xbt_long_nom);
xbt_time = double(xbt_time);
xbt_coord = double(xbt_coord); %input will be xbt_lat for zonal transects and xbt_long for meridional transects
xbt_depth = double(xbt_depth);
T = double(T);
S = double(S);

%Other variables
g=9.81;

%Bathymetry
topo_long = ncread('etopo6.nc','LON6');
topo_lat = ncread('etopo6.nc','LAT6101_1700');
topo_bath = ncread('etopo6.nc','BATH6');

%% Read in Argo climatology
%Load Argo values from high resolution grid
argo_mean2018_lat = double(ncread('RG_ArgoClim_33pfit_2019_mean.nc','LATITUDE'));
argo_mean2018_long = double(ncread('RG_ArgoClim_33pfit_2019_mean.nc','LONGITUDE'));
argo_mean2018_p = double(ncread('RG_ArgoClim_33pfit_2019_mean.nc','PRESSURE'));
argo_mean2018_T = ncread('RG_ArgoClim_33pfit_2019_mean.nc','ARGO_TEMPERATURE_MEAN');
argo_mean2018_S = ncread('RG_ArgoClim_33pfit_2019_mean.nc','ARGO_SALINITY_MEAN');

%Restrict Argo data range to region of interest
argo_long = argo_mean2018_long(argo_mean2018_long > xx(1) & argo_mean2018_long < xx(2));
argo_lat = argo_mean2018_lat(argo_mean2018_lat > yy(1) & argo_mean2018_lat < yy(2));
argo_T = argo_mean2018_T(argo_mean2018_long > xx(1) & argo_mean2018_long < xx(2),...
    argo_mean2018_lat > yy(1) & argo_mean2018_lat < yy(2),...
    :);
argo_S = argo_mean2018_S(argo_mean2018_long > xx(1) & argo_mean2018_long < xx(2),...
    argo_mean2018_lat > yy(1) & argo_mean2018_lat < yy(2),...
    :);
%permute Argo T and S so that it is lat x long x p
argo_T = permute(argo_T,[2 1 3]);
argo_S = permute(argo_S,[2 1 3]);
argo_P = argo_mean2018_p;

%Repeat coastal Argo values
if range(xbt_long_nom) > range(xbt_lat_nom) %zonal transect
    %Fill Argo NaN values at edges with the non-NaN T-S value one longitude grid point over
    argo_T_temp = argo_T;
    argo_S_temp = argo_S;
    for i=1:length(argo_lat)
        idx = find(isnan(argo_T(i,:,1))); %find NaN indices
        %for each NaN grid point
        for j=1:length(idx)
            %if the grid point to the E is non-NaN infill with these values
            if idx(j)~=length(argo_long) & ~isnan(argo_T(i,idx(j)+1,:))
                argo_T_temp(i,idx(j),:) = argo_T(i,idx(j)+1,:);
                argo_S_temp(i,idx(j),:) = argo_S(i,idx(j)+1,:);
                %else if the grid point to the W is non-NaN infill with these values
            elseif idx(j)~=1 & ~isnan(argo_T(i,idx(j)-1,:))
                argo_T_temp(i,idx(j),:) = argo_T(i,idx(j)-1,:);
                argo_S_temp(i,idx(j),:) = argo_S(i,idx(j)-1,:);
            end
        end
    end
    
elseif range(xbt_lat_nom) > range(xbt_long_nom) %meridional transect
    %Fill Argo NaN values at edges with the non-NaN T-S value one latitude grid point over
    argo_T_temp = argo_T;
    argo_S_temp = argo_S;
    for i=1:length(argo_long)
        idx = find(isnan(argo_T(:,i,1))); %find NaN indices
        %for each NaN grid point
        for j=1:length(idx)
            %if the grid point to the N is non-NaN infill with these values
            if idx(j)~=length(argo_lat) & ~isnan(argo_T(idx(j)+1,i,:))
                argo_T_temp(idx(j),i,:) = argo_T(idx(j)+1,i,:);
                argo_S_temp(idx(j),i,:) = argo_S(idx(j)+1,i,:);
                %else if the grid point to the S is non-NaN infill with these values
            elseif idx(j)~=1 & ~isnan(argo_T(idx(j)-1,i,:))
                argo_T_temp(idx(j),i,:) = argo_T(idx(j)-1,i,:);
                argo_S_temp(idx(j),i,:) = argo_S(idx(j)-1,i,:);
            end
        end
    end
end
%save values
argo_T = argo_T_temp;
argo_S = argo_S_temp;

%% -- Argo transect path correction --
%Use Argo climatology to correct for differences in the path of each HR-XBT
%transect compared to the nominal transect

%Interpolate Argo to HR-XBT depths
[Xi,Yi,Zi] = meshgrid(argo_long,argo_lat,xbt_depth); %create meshgrid for interpolation
argo_T_zint = interp3(argo_long,argo_lat,argo_P,argo_T,Xi,Yi,Zi);
argo_S_zint = interp3(argo_long,argo_lat,argo_P,argo_S,Xi,Yi,Zi);
%shallowest Argo data is at 2.5 dbar but shallowest xbt depth is 0 m so
%assume 2.5 dbar is representative of the surface
argo_T_zint(:,:,1) = argo_T(:,:,1);
argo_S_zint(:,:,1) = argo_S(:,:,1);


%Interpolate Argo to nominal transect
argo_T_nom = NaN(length(xbt_long_nom),length(xbt_depth));
argo_S_nom = NaN(length(xbt_long_nom),length(xbt_depth));
for z=1:length(xbt_depth)
    argo_T_nom(:,z) = interp2(argo_long,argo_lat,argo_T_zint(:,:,z),xbt_long_nom,xbt_lat_nom);
    argo_S_nom(:,z) = interp2(argo_long,argo_lat,argo_S_zint(:,:,z),xbt_long_nom,xbt_lat_nom);
end
%interpolate along nominal transect to infill any gaps
for z=1:length(xbt_depth)
    argo_T_nom(:,z) = interp1(xbt_long_nom(~isnan(argo_T_nom(:,z))),argo_T_nom(~isnan(argo_T_nom(:,z)),z),xbt_long_nom);
    argo_S_nom(:,z) = interp1(xbt_long_nom(~isnan(argo_S_nom(:,z))),argo_S_nom(~isnan(argo_S_nom(:,z)),z),xbt_long_nom);
end


%Adjust HR-XBT measurements onto nominal transect
T_xa_nom = NaN*T; %initialise
S_xa_nom = NaN*S;
if range(xbt_long_nom) > range(xbt_lat_nom) %Zonal transect
    for t=1:length(xbt_time)
        %interpolate Argo to HR-XBT transect
        argo_T_transect = NaN(length(xbt_long_nom),length(xbt_depth));
        argo_S_transect = NaN(length(xbt_long_nom),length(xbt_depth));
        for z=1:length(xbt_depth)
            argo_T_transect(:,z) = interp2(argo_long,argo_lat,argo_T_zint(:,:,z),xbt_long_nom,xbt_coord(:,t));
            argo_S_transect(:,z) = interp2(argo_long,argo_lat,argo_S_zint(:,:,z),xbt_long_nom,xbt_coord(:,t));
            %linearly interpolate to infill any gaps in Argo along the
            %transect, and extrapolate to include points with xbt data that
            %are missed due to the triangulation used by interp2
            idx = find(~isnan(argo_T_transect(:,z)));
            argo_T_transect(:,z) = interp1(xbt_long_nom(idx),argo_T_transect(idx,z),xbt_long_nom,'linear','extrap');
            argo_S_transect(:,z) = interp1(xbt_long_nom(idx),argo_S_transect(idx,z),xbt_long_nom,'linear','extrap');
        end
        
        %difference between Argo on HR-XBT transect and Argo on nominal transect
        T_diff = argo_T_transect - argo_T_nom;
        S_diff = argo_S_transect - argo_S_nom;
        
        %correct HR-XBT measurements onto nominal transect using Argo differencing
        T_xa_nom(:,:,t) = T(:,:,t) - T_diff;
        S_xa_nom(:,:,t) = S(:,:,t) - S_diff;
    end
    
elseif range(xbt_lat_nom) > range(xbt_long_nom) %Meridional transect
    for t=1:length(xbt_time)
        %interpolate Argo to HR-XBT transect
        argo_T_transect = NaN(length(xbt_long_nom),length(xbt_depth));
        argo_S_transect = NaN(length(xbt_long_nom),length(xbt_depth));
        for z=1:length(xbt_depth)
            argo_T_transect(:,z) = interp2(argo_long,argo_lat,argo_T_zint(:,:,z),xbt_coord(:,t),xbt_lat_nom);
            argo_S_transect(:,z) = interp2(argo_long,argo_lat,argo_S_zint(:,:,z),xbt_coord(:,t),xbt_lat_nom);
            %linearly interpolate to infill any gaps in Argo along the transect
            idx = find(~isnan(argo_T_transect(:,z)));
            argo_T_transect(:,z) = interp1(xbt_lat_nom(idx),argo_T_transect(idx,z),xbt_lat_nom);
            argo_S_transect(:,z) = interp1(xbt_lat_nom(idx),argo_S_transect(idx,z),xbt_lat_nom);
        end
        
        %difference between Argo on HR-XBT transect and Argo on nominal transect
        T_diff = argo_T_transect - argo_T_nom;
        S_diff = argo_S_transect - argo_S_nom;
        
        %correct HR-XBT measurements onto nominal transect using Argo differencing
        T_xa_nom(:,:,t) = T(:,:,t) - T_diff;
        S_xa_nom(:,:,t) = S(:,:,t) - S_diff;
    end
end

%% -- Find start and end longitudes --
%Take the edges of the nominal transect to be where at least 75% of data is
%still available (i.e. at most 25% of transects start after the starting
%point, and at most 25% of transects end before the ending point)
%nb. have only applied to zonal transects

%Start point
starting_coords = NaN*xbt_time;
for i=1:length(xbt_time)
    starting_coords(i) = xbt_long_nom(find(~all(isnan(T_xa_nom(:,:,i)')),1,'first'));
end
[A,B] = histcounts(starting_coords,xbt_long_nom-0.05,'Normalization','cdf');
start_coord = B(find(A>=0.75,1,'first'));
%End point
ending_coords = NaN*xbt_time;
for i=1:length(xbt_time)
    ending_coords(i) = xbt_long_nom(find(~all(isnan(T_xa_nom(:,:,i)')),1,'last'));
end
[C,D] = histcounts(ending_coords,xbt_long_nom+0.05,'Normalization','cdf');
end_coord = D(find(C>0.25,1,'first'))+0.1;

disp([start_coord+0.05 end_coord-0.05]) %display

%restrict lat, long, T, S to this new range (ensure double format)
T_nom = double(T_xa_nom(xbt_long_nom >= start_coord & xbt_long_nom <= end_coord,:,:));
S_nom = double(S_xa_nom(xbt_long_nom >= start_coord & xbt_long_nom <= end_coord,:,:));
long_nom = double(xbt_long_nom(xbt_long_nom >= start_coord & xbt_long_nom <= end_coord));
lat_nom = double(xbt_lat_nom(xbt_long_nom >= start_coord & xbt_long_nom <= end_coord));

%% -- Extend SH in depth using Argo  --
%Use Argo profiles to find linear relationship between temperature at a
%reference depth and SH relative to 1975 dbar
%nb. have only applied to zonal transects

% SH(z/1975) - SH(z/z_ref) ~ m(z_ref)*T_800 + c(z_ref), 0<=z<=z_ref
% SH(z/1975) ~ m(z)*T_800 + c(z), z_ref<z<=1975

%Read in Argo profiles
argo_prof_long = ncread(argo_fname,'LONGITUDE');
argo_prof_lat = ncread(argo_fname,'LATITUDE');
argo_prof_time = ncread(argo_fname,'TIME'); %days since 1986-01-01 00:00:00
argo_prof_p = ncread(argo_fname,'PRES');
argo_prof_T = ncread(argo_fname,'TEMP');
argo_prof_S = ncread(argo_fname,'PSAL');

argo_prof_time = double(datenum('01-Jan-1986') + argo_prof_time);

%identify +/-1.5 lat region around nominal transect
kink = find(long_nom>kink_long-0.05 & long_nom<kink_long+0.05);
long_nom_box = [long_nom(1), long_nom(kink), long_nom(end),...
    long_nom(end), long_nom(kink), long_nom(1),...
    long_nom(1)];
lat_nom_box = [lat_nom(1)+1.5, lat_nom(kink)+1.5, lat_nom(end)+1.5,...
    lat_nom(end)-1.5, lat_nom(kink)-1.5, lat_nom(1)-1.5,...
    lat_nom(1)+1.5];

%find profiles in region of interest
reg_idx = inpolygon(argo_prof_long,argo_prof_lat,long_nom_box,lat_nom_box);
argo_prof_long_reg = argo_prof_long(reg_idx);
argo_prof_lat_reg = argo_prof_lat(reg_idx);
argo_prof_time_reg = argo_prof_time(reg_idx);
argo_prof_p_reg = argo_prof_p(:,reg_idx);
argo_prof_T_reg = argo_prof_T(:,reg_idx);
argo_prof_S_reg = argo_prof_S(:,reg_idx);

%consider only profiles that extend to at least 1975 dbar
idx1975 = find(argo_prof_p(:,1) == 1975); %depth idx of the 1975 depth level
argo_prof_p_reg(isnan(argo_prof_T_reg)) = NaN; %for the Argo depth levels for which T is NaN, set P to NaN
z_idx = find(~isnan(argo_prof_p_reg(idx1975,:))); %find idx of profiles for which there is data at 1975 dbar
argo_prof_long_z = argo_prof_long_reg(z_idx);
argo_prof_lat_z = argo_prof_lat_reg(z_idx);
argo_prof_time_z = argo_prof_time_reg(z_idx);
argo_prof_p_z = argo_prof_p_reg(:,z_idx);
argo_prof_T_z = argo_prof_T_reg(:,z_idx);
argo_prof_S_z = argo_prof_S_reg(:,z_idx);

figure()
hold on
plot(argo_prof_long,argo_prof_lat,'.','MarkerSize',1,'DisplayName','Argo profiles')
plot(argo_prof_long_reg,argo_prof_lat_reg,'r.','MarkerSize',1,'DisplayName','Argo profiles (region)')
plot(argo_prof_long_z,argo_prof_lat_z,'g.','MarkerSize',1,'DisplayName','Argo profiles (1975-m)')
plot(long_nom,lat_nom,'k--','LineWidth',2,'DisplayName','Nominal transect')
plot(long_nom_box,lat_nom_box,'k','LineWidth',2,'DisplayName','Model region')
contour(topo_long,topo_lat,topo_bath',[-1000 -1000],'LineWidth',1,'LineColor',rgb('silver'),'DisplayName','1000-m isobath')
xlim(xx)
ylim(yy)
daspect([1 1 1])
xlabel('Longitude [\circE]')
ylabel('Latitude [\circN]')
box on
legend('Location','West')


%Compute Argo profile SH
%remove lowest level (so 1975 is the bottom)
argo_prof_p_1975 = argo_prof_p_z(1:end-1,:);
argo_prof_T_1975 = argo_prof_T_z(1:end-1,:);
argo_prof_S_1975 = argo_prof_S_z(1:end-1,:);

argo_p_1975 = argo_prof_p(1:end-1,1); %Argo depth level vector

%convert practical salinity to absolute salinity
SA = gsw_SA_from_SP(argo_prof_S_1975,argo_prof_p_1975,argo_prof_long_z,argo_prof_lat_z);
%convert in-situ temperature to conservative temperature
CT = gsw_CT_from_t(SA,argo_prof_T_1975,argo_prof_p_1975);
%calculate dynamic height relative to bottom p level (1975 dbar)
DH=NaN*CT; %[depth] x [# profiles]
for i=1:length(argo_prof_long_z)
    DH(:,i) = gsw_geo_strf_dyn_height(SA(:,i),CT(:,i),argo_prof_p_1975(:,i),argo_prof_p_1975(end,i));
end
%convert DH to SH
argo_prof_SH = DH/g;


%Compute Argo profile correlation between SH and reference temperature
ref_idx = find(argo_p_1975 == ref_depth); %find what depth idx is the reference depth level
argo_prof_Tref = argo_prof_T_1975(ref_idx,:); %argo profiles T at reference depth
argo_prof_SHref = argo_prof_SH - argo_prof_SH(ref_idx,:); %compute Argo profile SH relative to reference depth

r_SH=NaN*argo_p_1975;
for z=1:length(argo_p_1975)
    nan_idx = find(~isnan(argo_prof_SH(z,:)) & ~isnan(argo_prof_Tref)); %exclude NaNs
    if length(nan_idx) < length(argo_prof_SH)*0.5
        r = [NaN NaN]; %don't compute correlation if over half the profiles are NaN at a given depth
    elseif argo_p_1975(z) <= ref_depth
        r = corrcoef(argo_prof_SH(z,nan_idx)-argo_prof_SHref(z,nan_idx),argo_prof_Tref(nan_idx));
    else
        r = corrcoef(argo_prof_SH(z,nan_idx),argo_prof_Tref(nan_idx));
    end
    r_SH(z) = r(2);
end
figure()
plot(r_SH,-argo_p_1975,'LineWidth',2)
xlim([0 1])
title('R - reference temperature and SH')
box on
text(0.3,-200,num2str(nanmean(r_SH)))

%Build least squares model from profiles
model=NaN(length(argo_p_1975),2);
for z=1:length(argo_p_1975)
    nan_idx = find(~isnan(argo_prof_Tref)); %exclude NaNs while building model
    G = [ones(length(argo_prof_Tref(nan_idx)),1), argo_prof_Tref(nan_idx)'];
    if argo_p_1975(z) <= ref_depth
        d=argo_prof_SH(z,nan_idx)' - argo_prof_SHref(z,nan_idx)';
    else
        d=argo_prof_SH(z,nan_idx)';
    end
    model(z,:) = G\d; %[m(1) is the intercept; m(2) is the gradient]
end
%regression coeffecients at reference depth are used for all depths shallower
model(argo_p_1975 < ref_depth,:) = repmat(model(ref_idx,:),length(find(argo_p_1975 < ref_depth)),1);


%Apply model to XBT transects
SH_nom_1975 = NaN(length(argo_p_1975),length(long_nom),length(xbt_time)); %[depth] x [long] x [time]
for t=1:length(xbt_time )
    %interpolate XBT depths to Argo profile depth levels
    T_use = NaN(length(long_nom),length(argo_p_1975));
    S_use = T_use;
    for i=1:length(long_nom)
        T_use(i,:) = interp1(xbt_depth,T_nom(i,:,t),argo_p_1975);
        S_use(i,:) = interp1(xbt_depth,S_nom(i,:,t),argo_p_1975);
    end
    T_xref = T_use(:,argo_p_1975 == ref_depth); %T at reference depth
   
    %convert practical salinity to absolute salinity
    SA = gsw_SA_from_SP(S_use',argo_p_1975,long_nom,lat_nom);
    %convert in-situ temperature to conservative temperature
    CT = gsw_CT_from_t(SA,T_use',argo_p_1975);
    %calculate dynamic height relative to reference depth
    DH=NaN*CT;
    for i=1:length(long_nom)
        DH(:,i) = gsw_geo_strf_dyn_height(SA(:,i),CT(:,i),argo_p_1975,argo_p_1975(argo_p_1975 == ref_depth));
    end
    %convert DH to SH
    SH_xref = DH/g;
    
    %Compute XBT SH relative to 1975 dbar using least squares model
    for z=1:length(argo_p_1975)
        if argo_p_1975(z) <= ref_depth
            SH_nom_1975(z,:,t) = model(z,1) + model(z,2)*T_xref + SH_xref(z,:)';
        else
            SH_nom_1975(z,:,t) = model(z,1) + model(z,2)*T_xref;
        end
    end
end


%Adjust SH to be relative to bathymetry or bottom p level (1975 dbar), whichever is shallower
%interpolate bathymetry onto nominal transect
topo_transect = -1*interp2(topo_lat,topo_long,topo_bath,lat_nom,long_nom);
%for each position along transect find which pressure level is closest to the bathymetry (assuming m is equivalent to dbar)
bath_idx = NaN*long_nom;
SH_nom = SH_nom_1975;
for i=1:length(long_nom)
    [~,bath_idx(i)] = min(abs(argo_p_1975-topo_transect(i)));
    SH_nom(:,i,:) = SH_nom_1975(:,i,:) - SH_nom_1975(bath_idx(i),i,:);
end


%Compute correlation of Argo profile SH with estimated Argo profile SH to test model
r_profiles=NaN*argo_p_1975;
SH_new = NaN*argo_prof_SH;
for z=1:length(argo_p_1975)
    %apply least squares model
    if argo_p_1975(z) <= ref_depth
        SH_new(z,:) = model(z,1) + model(z,2)*argo_prof_Tref + argo_prof_SHref(z,:);
    else
        SH_new(z,:) = model(z,1) + model(z,2)*argo_prof_Tref;
    end
    %calculate correlation coefficient
    nan_idx = find(~isnan(argo_prof_SH(z,:)) & ~isnan(SH_new(z,:)));
    if length(nan_idx) < length(argo_prof_SH)*0.5 %don't compute correlation if over half the profiles are NaN at a given depth
        r = [NaN NaN];
    else
        r = corrcoef(SH_new(z,nan_idx),argo_prof_SH(z,nan_idx));
        r_profiles(z) = r(2);
    end
end
figure()
plot(r_profiles,-argo_p_1975,'LineWidth',2)
xlim([0 1])
title('R - full model profiles')
box on

%% -- Use satellite altimetry to increase temporal resolution --
%Read in SLA
sla = ncread(sla_fname,'sla'); %referenced to 1993-2012 period %[long x lat x time]
time_sat = ncread(sla_fname,'time');
lat_sat = ncread(sla_fname,'latitude');
long_sat = ncread(sla_fname,'longitude');

%time is days since 1950-01-01
time_sat = double(time_sat + datenum('01-Jan-1950'));

%Interpolate SLA to nominal transect
sla_nom = NaN(length(long_nom),length(time_sat));
for i=1:length(time_sat)
    sla_nom(:,i) = interp2(lat_sat,long_sat,sla(:,:,i),lat_nom,long_nom);
end

%find SLA temporal anomaly over the HR-XBT transect time period (typically 2004-2019)
xbt_period_idx = find(time_sat >= datenum(dateshift(datetime(xbt_time(1),'ConvertFrom','Datenum'),'start','year')) &...
    time_sat <= datenum(dateshift(datetime(xbt_time(end),'ConvertFrom','Datenum'),'end','year')));

%compute SLA temporal anomaly
sla_nom_a = sla_nom - mean(sla_nom(:,xbt_period_idx),2);

%remove SLA trend over the HR-XBT transect time period
sla_nom_a_orig = sla_nom_a;
for i=1:length(long_nom)
    %compute linear fit coefficients over HR-XBT period
    sla_coef = polyfit(time_sat(xbt_period_idx),sla_nom_a(i,xbt_period_idx),1);
    %apply linear fit coefficients over full 2004-2019 period
    sla_trend = polyval(sla_coef,time_sat);
    %subtract trend from SLA
    sla_nom_a(i,:) = sla_nom_a_orig(i,:) - sla_trend';
end

%Read in xbt profile .dat files
%specify XBT .dat files
dat_path_name = ['XBT\dat\',trans_num];

%read in all .dat files
dotdat = dir(strcat(dat_path_name,'\*.dat'));
%initiate arrays (padded with excess NaNs) %[file x length]
xbt_prof_date = NaT(length(dotdat),400);
xbt_prof_lat = NaN(length(dotdat),400);
xbt_prof_long = NaN(length(dotdat),400);
%read each file
for i = 1:length(dotdat)
    try
        A = readtable(dotdat(i).name,'Format','%f %f %f %{dd/MM/yy}D %{hh:mm:ss}T %f %f %f %f %f');
    catch
        A = readtable(dotdat(i).name,'Format','%f %f %f %{dd/MM/yy}D %{hh:mm:ss}T %f %f %f %f %f %s');
    end
    date1 = A{:,4};
    lat1 = A{:,6};
    long1 = A{:,7};
    drop1 = A{:,9};
    %keep only good drops ("1 or 2 for good drop")
    bad_drops = find(drop1<1 | drop1>2);
    date1(bad_drops) = NaT;
    lat1(bad_drops) = NaN;
    long1(bad_drops) = NaN;
    drop1(bad_drops) = NaN;
    %correct dates if time is not read as being in 2000's
    if year(min(date1)) < 1000
        date1 = date1+years(2000);
    end
    if isa(date1,'datetime') == 0
        'Error: Time input read incorrectly'
        return
    end
    %store in array
    xbt_prof_date(i,1:length(date1)) = date1;
    xbt_prof_lat(i,1:length(lat1)) = lat1;
    xbt_prof_long(i,1:length(long1)) = long1;
end


%Find date of xbt profiles closest to the region of maximum sla variance
%find longitude of maximum SLA variance
sla_nom_var = var(sla_nom,0,2);
[~,I] = max(sla_nom_var);
max_sla_var_long = long_nom(I);

[~,I] = min(abs(xbt_prof_long-max_sla_var_long),[],2); %find index of closest longitude for each occupation
xbt_prof_long_maxvar = NaN*I;
xbt_prof_date_maxvar = NaT(size(I));
for i=1:length(I)
    xbt_prof_long_maxvar(i) = xbt_prof_long(i,I(i));
    xbt_prof_date_maxvar(i) = xbt_prof_date(i,I(i));
end

%Restrict to just the 2004-2019 period
date_idx = find(xbt_prof_date_maxvar>=datetime('01-Jan-2004') & xbt_prof_date_maxvar<=datetime('31-Dec-2019'));
xbt_prof_date_maxvar_04 = xbt_prof_date_maxvar(date_idx);

%Remove dates coninciding with removed transects
%keep the xbt profile dates if there is a corresponding date from the xbt transect timestamp within +/-30 days)
xbt_dates = xbt_prof_date_maxvar_04;
for t=1:length(xbt_prof_date_maxvar_04)
    if isempty(find(xbt_time > datenum(xbt_prof_date_maxvar_04(t))-30 & xbt_time < datenum(xbt_prof_date_maxvar_04(t))+30))
        xbt_dates(t) = NaT;
    end
end
xbt_dates = sort(datenum(rmmissing(xbt_dates)));


%Match SLA date with xbt profile date
time_sat_match_idx = NaN*xbt_dates;
count=1;
for t=1:length(time_sat)
    if ismember(time_sat(t),xbt_dates)
        time_sat_match_idx(count)=t;
        count=count+1;
    end
end
time_sat_x = time_sat(time_sat_match_idx);

%Average SLA over 7 days centred on the XBT profile date
sla_nom_a_x = NaN(length(long_nom),length(xbt_time));
for i=1:length(long_nom)
    for t=1:length(time_sat_match_idx)
        sla_nom_a_x(i,t) = mean(sla_nom_a(i,time_sat_match_idx(t)-3:time_sat_match_idx(t)+3));
    end
end


%Linear regression of SH' at the surface with SLA'
%Find XBT+Argo SH temporal anomaly
SH_nom_mean = nanmean(SH_nom,3);
SH_nom_a = SH_nom - SH_nom_mean; %[depth x long x time]

%remove SH trend from HR-XBT+Argo SH
SH_nom_a_orig = SH_nom_a;
%trend coefficients will be applied to 2004-2019 monthly time series to be added back in
SH_nom_trend_coeff = NaN(length(argo_p_1975),length(long_nom),2); %[depth x long x coefficients]
for i=1:length(long_nom)
    for z=1:length(argo_p_1975)
        %find NaNs
        idx = ~isnan(squeeze(SH_nom_a_orig(z,i,:)));
        %compute linear fit coefficients over HR-XBT period and save coefficients
        p = polyfit(xbt_dates(idx),SH_nom_a_orig(z,i,idx),1);
        SH_nom_trend_coeff(z,i,:) = p;  %[:,:,1 is the gradient; :,:,2 is the intercept]
        %apply linear fit coefficients over HR-XBT period to compute trend
        SH_trend_x = polyval(p,xbt_dates);
        %remove trend from de-meaned SH
        SH_nom_a(z,i,:) = squeeze(SH_nom_a_orig(z,i,:)) - SH_trend_x;
    end
end

%vector of the surface values
SHa_surf = squeeze(SH_nom_a(1,:,:));
SHa_vec = SHa_surf(:);
SLAa_vec = sla_nom_a_x(:);
surf_vals = ~isnan(SHa_vec); %indices of non-NaN values

%Linear regression
[p] = polyfit(SLAa_vec(surf_vals),SHa_vec(surf_vals),1); %[p(1) is the gradient; p(2) is the intercept]
regression_vals = polyval(p,SLAa_vec);

r = corrcoef(regression_vals(surf_vals),SHa_vec(surf_vals)); %correlation coefficient

%Plot
figure()
plot(SLAa_vec,SHa_vec,'b.')
hold on
plot(SLAa_vec,regression_vals,'r.')
xlim([-1 1])
ylim([-1 1])
xlabel('SLA temporal anomaly [m]')
ylabel('Surface SH temporal anomaly [m]')
text(-0.7,0.7,{['r = ', num2str(round(r(2),3))],['SH = ', num2str(round(p(1),3)), ' * SLA + ', num2str(round(p(2),3))]})


%Build regression model using 7-day averaged daily SLA'
sla_model = NaN(length(argo_p_1975),length(long_nom),2); %[depth x long x coefficients]
for i=1:length(long_nom)
    for z=1:length(argo_p_1975)
        idx = ~isnan(squeeze(SH_nom_a(z,i,:))); %find NaNs
        
        %linear regression
        sla_model(z,i,:) = polyfit(sla_nom_a_x(i,idx)',squeeze(SH_nom_a(z,i,idx)),1); %[:,:,1 is the gradient; :,:,2 is the intercept]
    end
end


%Monthly average daily SLA 
%initialise
time_monthly = datenum(2004,1:192,15);
sla_nom_monthly_a = NaN(length(long_nom),length(time_monthly));

%monthly average daily data
[y_monthly,m_monthly,~] = ymd(datetime(time_monthly,'ConvertFrom','datenum'));
[y_daily,m_daily,~] = ymd(datetime(time_sat,'ConvertFrom','datenum'));
for t=1:length(time_monthly)
    sla_nom_monthly_a(:,t) = mean(sla_nom_a(:,(y_daily == y_monthly(t) & m_daily == m_monthly(t))),2);   
end


%Apply regression model to monthly SLA'
SH_nom_monthly_a = NaN(length(argo_p_1975),length(long_nom),length(time_monthly)); %[depth x long x time]
%apply regression at each point at each depth
for i=1:length(long_nom)
    for z=1:length(argo_p_1975)
        SH_nom_monthly_a(z,i,:) = sla_nom_monthly_a(i,:)*sla_model(z,i,1) + sla_model(z,i,2);
    end
end


%Add HR-XBT+Argo SH trend back in
%monthly time series index for range of years which have HR-XBT measurements  
trend_idx = find(time_monthly >= datenum(dateshift(datetime(xbt_dates(1),'ConvertFrom','Datenum'),'start','year')) &...
    time_monthly <= datenum(dateshift(datetime(xbt_dates(end),'ConvertFrom','Datenum'),'end','year')));
%add HR-XBT+Argo SH trend back in
SH_nom_monthly_a_wtrend = NaN*SH_nom_monthly_a; %initialise
for i=1:length(long_nom)
    for z=1:length(argo_p_1975)
        %set trend as 0 outside of HR-XBT sampling period
        SH_trend = zeros(length(time_monthly),1);
        %compute linear trend from linear fit coefficients only over HR-XBT period
        SH_trend(trend_idx) = polyval(squeeze(SH_nom_trend_coeff(z,i,:)),time_monthly(trend_idx));
        %add trend to monthly SH anomaly 
        SH_nom_monthly_a_wtrend(z,i,:) = squeeze(SH_nom_monthly_a(z,i,:)) + SH_trend;
    end
end
%Add mean SH back in
SH_nom_monthly = SH_nom_monthly_a_wtrend + SH_nom_mean;


%Linear regression of monthly SH' at the surface and monthly SLA'
%vector of all the surface values
SHa_surf = squeeze(SH_nom_monthly_a(1,:,:));
SHa_vec = SHa_surf(:);
SLAa_vec = sla_nom_monthly_a(:);

%Linear regression
[p] = polyfit(SLAa_vec,SHa_vec,1); %[p(1) is the gradient; p(2) is the intercept]
regression_vals = polyval(p,SLAa_vec);

r = corrcoef(regression_vals,SHa_vec); %correlation coefficient

%Plot
figure()
plot(SLAa_vec,SHa_vec,'b.')
hold on
plot(SLAa_vec,regression_vals,'r.')
xlim([-1 1])
ylim([-1 1])
xlabel('Monthly SLA temporal anomaly [m]')
ylabel('Monthly surface SH temporal anomaly [m]')
text(-0.7,0.7,{['r = ', num2str(round(r(2),3))],['SH = ', num2str(round(p(1),3)), ' * SLA + ', num2str(round(p(2),3))]})

%Compute correlation between XBT+Argo SH' and XBT+Argo+Altimetry SH'
%match monthly sla and xbt months
[y_xbt,m_xbt] = ymd(datetime(xbt_dates,'ConvertFrom','datenum'));
[y_sat,m_sat] = ymd(datetime(time_monthly,'ConvertFrom','datenum'));
month_match_idx = NaN(length(xbt_time),1);
for t=1:length(month_match_idx)
    month_match_idx(t) = find(m_xbt(t) == m_sat & y_xbt(t) == y_sat);
end

r_section = NaN(length(argo_p_1975),length(long_nom));
%compute correlation
for i=1:length(long_nom)
    for z=1:length(argo_p_1975)
        idx = ~isnan(squeeze(SH_nom_a(z,i,:))); %find NaNs
        
        %compute correlation between SH from altimetry and XBT+Argo SH
        r = corrcoef(SH_nom_monthly_a(z,i,month_match_idx(idx)),SH_nom_a(z,i,idx));
        r_section(z,i) = r(2);
    end
end

figure()
contourf(long_nom,-argo_p_1975,abs(r_section),'LineColor','none')
colorbar
colormap(parula(10))
caxis([0 1])
ylabel(colorbar,'abs(r)')

%% -- Compute geostrophic velocity --
DH_nom = SH_nom_monthly*g; %convert SH to DH
gvel_nom = NaN*DH_nom;
for t=1:length(time_monthly)
    idx = find(all(~isnan(DH_nom(:,:,t)),1)); %only considering points with no-NaN
    gvel_nom(:,idx,t) = calc_gvel(DH_nom(:,idx,t),argo_p_1975,long_nom(idx),lat_nom(idx),bath_idx(idx));
end

end %End function