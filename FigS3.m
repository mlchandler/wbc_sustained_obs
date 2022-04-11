% Mitchell Chandler, SIO
% Last updated: 11/04/2022

%% Read in monthly PX30 data
time_monthly = ncread('px30_velocity.nc','time')+datenum('01-Jan-1970');
px30_lat_nom = ncread('px30_velocity.nc','latitude');
px30_long_nom = ncread('px30_velocity.nc','longitude');
argo_depth = ncread('px30_velocity.nc','depth');
px30_gvel_LKM = ncread('px30_velocity.nc','vel');
px30_gvel_LNM = ncread('px30_velocity.nc','gvel_LNM');

%% Read in mooring data from Sloyan et al. 2016
load for_nathalie_norot
timed = timed + datenum('01-Jan-1900');

%% --- LKM ---
%% Rotate u and v velocity to compute net transport over the full mooring array normal to transect PX30
%find angle of mooring array from horizontal - E=0, N=90, S=-90
[~,phaseangle_moor] = sw_dist([-27.34 -27.33],[x_hr(1) x_hr(end)]);

%only consider PX30 measurements over similar longitude range as mooring array
long_idx = find(px30_long_nom>=153.75 & px30_long_nom<=155.3);

%find angle of PX30 from horizontal - E=0, N=90, S=-90
[~,phaseangle_px30] = sw_dist([px30_lat_nom(long_idx(1)) px30_lat_nom(long_idx(end))],[px30_long_nom(long_idx(1)) px30_long_nom(long_idx(end))]);

theta = phaseangle_px30-phaseangle_moor; %px30 is rotated this much counterclockwise from the mooring array

%rotate velocity from normal to mooring array to normal to px30 
cross_transect_moor_vel = U_across*sind(theta) + U_along*cosd(theta); %y' =  x sin(theta) + y cos(theta) 

%restrict to velocity between 0--1980-m
depth_idx = find(depth <= 1980); 
depth2 = depth(depth_idx);
ct_moor_vel2 = cross_transect_moor_vel(:,depth_idx,:); %[t,z,x]
ct_moor_vel2 = permute(ct_moor_vel2,[2 3 1]); %permute to be [z,x,t]

%depth-integrated velocity
[ct_moor_DIV] = calc_depth_int_vel(x_hr,depth2,ct_moor_vel2,timed);

%transport
lat_approx = interp1(x_hr([1 end]),[-27.34 -27.10],x_hr); %estimate latitude using first and last values from Sloyan et al. 2016
cumtransport_moor = calc_cumulative_transport(x_hr',lat_approx',ct_moor_DIV,timed);
ct_moor_transport2 = cumtransport_moor(end,:);

%monthly-average
% [datestr(timed(1)); datestr(timed(end))]
timed_monthly = datenum(2012,4:20,15);
tsouth_monthly = NaN*timed_monthly;

[y_monthly,m_monthly,~] = ymd(datetime(timed_monthly,'ConvertFrom','datenum'));
[y_moor,m_moor,~] = ymd(datetime(timed,'ConvertFrom','datenum'));
ct_transport_monthly = NaN*timed_monthly;
for t=1:length(timed_monthly)
    ct_transport_monthly(t) = mean(ct_moor_transport2(y_moor == y_monthly(t) & m_moor == m_monthly(t)));
end

%% Compute XAA transport over mooring array longitude
%restrict longitude range
px30_new_vel = px30_gvel_LKM(:,long_idx,:);
px30_new_lat = px30_lat_nom(long_idx);
px30_new_long = px30_long_nom(long_idx);

%depth-integrated velocity
[monthly_div] = calc_depth_int_vel(px30_new_long,argo_depth,px30_new_vel,time_monthly);

%transport
cumtransport_px30 = calc_cumulative_transport(px30_new_long,px30_new_lat,monthly_div,time_monthly);
px30_transport = cumtransport_px30(end,:);

%% 3-month running mean for absolute velocities
%(same as used for producing annual cycle)
px30_transport_smooth = movmean(px30_transport,3,'Endpoints','fill');
moored_transport_smooth = movmean(ct_transport_monthly,3,'Endpoints','fill');

%% --- LNM ---
%% Re-reference cross-transect mooring velocity to LNM
%find last non-nan value in each column
idx = zeros(length(depth2),length(x_hr));
mean_ct_moor_vel2 = nanmean(ct_moor_vel2,3);
for i=1:length(x_hr)
    V = mean_ct_moor_vel2(:,i);
    A = find(~isnan(V),1,'last');
    idx(A,i) = 1;
end
IDX = logical(idx);

%re-reference
ct_moor_vel2_LNM = NaN*ct_moor_vel2;
for t=1:length(timed)
    V = ct_moor_vel2(:,:,t);
    ct_moor_vel2_LNM(:,:,t) = V - V(IDX)';
end

%% Compute cross-transect transport from LNM
%depth-integrated velocity
[ct_moor_DIV_LNM] = calc_depth_int_vel(x_hr,depth2,ct_moor_vel2_LNM,timed);
%transport
cumtransport_moor = calc_cumulative_transport(x_hr',lat_approx',ct_moor_DIV_LNM,timed);
ct_moor_transport2_LNM = cumtransport_moor(end,:);

%monthly-average
ct_transport_monthly_LNM = NaN*timed_monthly;
for t=1:length(timed_monthly)
    ct_transport_monthly_LNM(t) = mean(ct_moor_transport2_LNM(y_moor == y_monthly(t) & m_moor == m_monthly(t)));
end

%% XAA over mooring array
%restrict longitude range
px30_new_vel_LNM = px30_gvel_LNM(:,long_idx,:);
%depth-integrated velocity
[monthly_div_LNM] = calc_depth_int_vel(px30_new_long,argo_depth,px30_new_vel_LNM,time_monthly);
%transport
cumtransport_px30 = calc_cumulative_transport(px30_new_long,px30_new_lat,monthly_div_LNM,time_monthly);
px30_transport_LNM = cumtransport_px30(end,:);

%% --- Plot ---
fsize = 11;

sc = 7E-2;

figure('color','w')
clf

%LNM
subplot(2,1,1)
hold on
%transport:
plot(timed,ct_moor_transport2_LNM,'LineWidth',1,'DisplayName','IMOS (daily)')
plot(timed_monthly,ct_transport_monthly_LNM,'Color',[0 0.4470 0.7410],'LineWidth',2.5,'DisplayName','IMOS (monthly)')
plot(time_monthly,px30_transport_LNM,'Color',[0.8500 0.3250 0.0980],'LineWidth',2.5,'DisplayName','XAA (monthly)')
%axis:
xlim([datenum('15-Jan-2012') datenum('15-Jan-2014')])
xticks(datenum(2012,1:4:25,15))
datetick('x','mmm yyyy','keeplimits','keepticks')
ylim([-60 15])
yticks(-60:15:15)
ylabel('Transport [Sv]')
%other:
legend('Location','S','Orientation','Horizontal','EdgeColor','w')
box on
YL = ylim;
ypos = max(YL) + range(YL)*sc;
text(datenum('15-Jan-2014'),ypos,'(a) LNM','FontSize',fsize,'HorizontalAlignment','Right')
ax = gca;
ax.FontSize = fsize;


%LKM
subplot(2,1,2)
hold on
%transport:
plot(timed,ct_moor_transport2,'LineWidth',1,'DisplayName','IMOS (daily)')
plot(timed_monthly,moored_transport_smooth,'Color',[0 0.4470 0.7410],'LineWidth',2.5,'DisplayName','IMOS (monthly)')
plot(time_monthly,px30_transport_smooth,'Color',[0.8500 0.3250 0.0980],'LineWidth',2.5,'DisplayName','XAA (monthly)')
%axis:
xlim([datenum('15-Jan-2012') datenum('15-Jan-2014')])
xticks(datenum(2012,1:4:25,15))
datetick('x','mmm yyyy','keeplimits','keepticks')
ylim([-60 15])
yticks(-60:15:15)
ylabel('Transport [Sv]')
%other:
legend('Location','S','Orientation','Horizontal','EdgeColor','w')
box on
YL = ylim;
ypos = max(YL) + range(YL)*sc;
text(datenum('15-Jan-2014'),ypos,'(b) LKM','FontSize',fsize,'HorizontalAlignment','Right')
ax = gca;
ax.FontSize = fsize;



