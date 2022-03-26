% Mitchell Chandler, SIO
% Last updated: 26/03/2022

%% Read in monthly PX40 velocity
px40_time = ncread('px40_velocity.nc','time')+datenum('01-Jan-1970');
px40_lat = ncread('px40_velocity.nc','latitude');
px40_long = ncread('px40_velocity.nc','longitude');
px40_depth = ncread('px40_velocity.nc','depth');
px40_gvel = ncread('px40_velocity.nc','vel');
px40_gvel_LNM = ncread('px40_velocity.nc','gvel_LNM');

%% Read in KESS
load OImapBottomPrsVn.mat %bottom reference values and OI lat-long grid
load mapday1475.mat %for getting matrix sizes

%% Mask reference values where error is large
%Mask out or exclude regions where the mapped bottom pressure error
%(errprs) exceeds 0.02 dbar.
errprs_cutoff = 0.02;

u_ref = OIprs.u;
u_ref(OIprs.errprs > errprs_cutoff) = NaN;
v_ref = OIprs.v;
v_ref(OIprs.errprs > errprs_cutoff) = NaN;

%% Read in all KESS gridded data
grid_file = dir('[redacted file name]');
%initialise variables
time = NaN(size(grid_file)); %dates of measurements
pressure = mapday.pressure; %pressure levels
velocity = NaN(length(px40_long),length(pressure),length(time)); %velocity normal to transect
for gf = 1:length(grid_file)
    load(grid_file(gf).name) %load file
    oidate_idx = find(mapday.dd == OIprs.dd); %find OI idx corresponding to this date
    if ~isempty(oidate_idx) %proceed only if date match between grid and reference
        time(gf) = mapday.dd + datenum('01-Jan-2004'); %save time, referenced to 1 January 2004
        
        %add the reference velocities to the gridded velocities
        u = mapday.u + u_ref(:,:,oidate_idx);
        v = mapday.v + v_ref(:,:,oidate_idx);
        
        %linearly interpolate absolute velocities onto nominal transect
        nom_u = NaN(length(px40_long),length(pressure));
        nom_v = nom_u;
        for z=1:length(pressure)
            nom_u(:,z) = interp2(OIprs.lon,OIprs.lat,u(:,:,z),px40_long,px40_lat);
            nom_v(:,z) = interp2(OIprs.lon,OIprs.lat,v(:,:,z),px40_long,px40_lat);
        end
        
        %linearly interpolate along-track to infill any gaps
        for z=1:length(pressure)
            ii = find(~isnan(nom_u(:,z)));
            nom_u(:,z) = interp1(px40_long(ii),nom_u(ii,z),px40_long);
            nom_v(:,z) = interp1(px40_long(ii),nom_v(ii,z),px40_long);
        end
        
        %rotate velocities to obtain the component normal to the transect
        [~,a] = sw_dist(px40_lat,px40_long); %angle of transect from E where N is +ve and S is -ve
        theta = mean(a);
        R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; %rotation matrix...
        %rotates points in the x-y plane counterclockwise through angle theta wrt the x-axis.
        N = [0; 1]; %northwards unit vector
        unit_vec = R*N; %unit vector orthogonal to the transect
        vv=NaN*nom_u;
        for i=1:length(px40_long)
            for z=1:length(pressure)
                vv(i,z) = dot(unit_vec,[nom_u(i,z); nom_v(i,z)]); %dot product to project velocity onto vector normal to transect
            end
        end
        %save velocity
        velocity(:,:,gf) = vv;
    end
end

%remove times that are NaN
velocity(:,:,isnan(time)) = [];
time(isnan(time)) = [];

%linearly interpolate to fill the one 1.5 day gap
time_full = [time(1):0.5:time(end)]';
v_full = NaN(length(px40_long),length(pressure),length(time_full));
for x=1:length(px40_long)
    for z=1:length(pressure)
    v_full(x,z,:) = interp1(time,squeeze(velocity(x,z,:)),time_full);
    end
end

%use mean to find where there is at least one NaN
mean_velocity = mean(v_full,3);
keep_idx = find(~isnan(mean_velocity(:,1)));
%keep values for which the mapped values are always good 
KESS_vel = v_full(keep_idx,:,:);
KESS_time = time_full;
KESS_p = pressure;
KESS_long_nom = px40_long(keep_idx); 
KESS_lat_nom = px40_lat(keep_idx);


%% Match PX40 and KESS regions
%Restrict px40 to same time and longitude range as KESS
%monthly:
m_time_idx = px40_time>=KESS_time(1) & px40_time<=KESS_time(end);
m_long_idx = px40_long>=KESS_long_nom(1) & px40_long<=KESS_long_nom(end);
monthly_gvel_LNM = px40_gvel_LNM(:,m_long_idx,m_time_idx);
monthly_gvel_LKM = px40_gvel(:,m_long_idx,m_time_idx);
monthly_time = px40_time(m_time_idx);
monthly_long = px40_long(m_long_idx);
monthly_lat = px40_lat(m_long_idx);

%Restrict KESS to same depth range as XAA (KESS upper 2000, XAA upper 1975)
k_pres_idx = KESS_p <= 2000;
kess_vel_LKM = KESS_vel(:,k_pres_idx,:);
kess_vel_LNM = kess_vel_LKM - kess_vel_LKM(:,end,:);
kess_p = KESS_p(k_pres_idx);
%and permute for [z,x,t] same as XAA
kess_vel_LKM = permute(kess_vel_LKM,[2 1 3]);
kess_vel_LNM = permute(kess_vel_LNM,[2 1 3]);

%% Compute transport
%LNM:
%Depth-integrated velocity
[monthly_DI_gvel_LNM] = calc_depth_int_vel(monthly_long,px40_depth,monthly_gvel_LNM,monthly_time);
[kess_DI_gvel_LNM] = calc_depth_int_vel(KESS_long_nom,kess_p,kess_vel_LNM,KESS_time);
%Transport
monthly_cumtransport_LNM = calc_cumulative_transport(monthly_long,monthly_lat,monthly_DI_gvel_LNM,monthly_time);
monthly_transport_LNM = monthly_cumtransport_LNM(end,:);
kess_cumtransport_LNM = calc_cumulative_transport(KESS_long_nom,KESS_lat_nom,kess_DI_gvel_LNM,KESS_time);
kess_transport_LNM = kess_cumtransport_LNM(end,:);

%LKM:
%Depth-integrated velocity
[monthly_DI_gvel_LKM] = calc_depth_int_vel(monthly_long,px40_depth,monthly_gvel_LKM,monthly_time);
[kess_DI_gvel_LKM] = calc_depth_int_vel(KESS_long_nom,kess_p,kess_vel_LKM,KESS_time);
%Transport
monthly_cumtransport_LKM = calc_cumulative_transport(monthly_long,monthly_lat,monthly_DI_gvel_LKM,monthly_time);
monthly_transport_LKM = monthly_cumtransport_LKM(end,:);
kess_cumtransport_LKM = calc_cumulative_transport(KESS_long_nom,KESS_lat_nom,kess_DI_gvel_LKM,KESS_time);
kess_transport_LKM = kess_cumtransport_LKM(end,:);

%% Monthly-average KESS transport
%initialise
kess_monthly_transport_LNM = NaN*monthly_transport_LNM;
kess_monthly_transport_LKM = NaN*monthly_transport_LKM;


[y_monthly,m_monthly,~] = ymd(datetime(monthly_time,'ConvertFrom','datenum'));
[y_k,m_k,~] = ymd(datetime(KESS_time,'ConvertFrom','datenum'));
for t=1:length(monthly_time)
    kess_monthly_transport_LNM(t) = mean(kess_transport_LNM(y_k == y_monthly(t) & m_k == m_monthly(t)));
    kess_monthly_transport_LKM(t) = mean(kess_transport_LKM(y_k == y_monthly(t) & m_k == m_monthly(t)));
end

%% Plot
fsize = 11;

figure()
clf

%LNM
subplot(2,1,1)
hold on
%transport:
plot(KESS_time,kess_transport_LNM,'LineWidth',1,'DisplayName','KESS (12-hourly)')
plot(monthly_time,kess_monthly_transport_LNM,'Color',[0 0.4470 0.7410],'LineWidth',2.5,'DisplayName','KESS (monthly)')
plot(monthly_time,monthly_transport_LNM,'Color',[0.8500 0.3250 0.0980],'LineWidth',2.5,'DisplayName','XAA (monthly)')
%axis:
ylim([-100 100])
ylabel('Transport [Sv]')
xlim([datenum('01-Jun-2004') datenum('01-Oct-2005')])
xticks(datenum(2004,6:4:6+16,1))
datetick('x','mmm yyyy','keeplimits','keepticks')
%other:
legend('Location','NW')
legend('boxoff')
box on
text(KESS_time(end),115,'(a) LNM','FontSize',fsize,'HorizontalAlignment','Right')
ax = gca;
ax.FontSize = fsize;

%LKM
subplot(2,1,2)
hold on
%transport:
plot(KESS_time,kess_transport_LKM,'LineWidth',1,'DisplayName','KESS (12-hourly)')
plot(monthly_time,kess_monthly_transport_LKM,'Color',[0 0.4470 0.7410],'LineWidth',2.5,'DisplayName','KESS (monthly)')
plot(monthly_time,monthly_transport_LKM,'Color',[0.8500 0.3250 0.0980],'LineWidth',2,'DisplayName','XAA (monthly)')
%axis:
ylim([-100 100])
ylabel('Transport [Sv]')
xlim([datenum('01-Jun-2004') datenum('01-Oct-2005')])
xticks(datenum(2004,6:4:6+16,1))
datetick('x','mmm yyyy','keeplimits','keepticks')
%other:
legend('Location','NW')
legend('boxoff')
box on
text(KESS_time(end),115,'(b) LKM','FontSize',fsize,'HorizontalAlignment','Right')
ax = gca;
ax.FontSize = fsize;


