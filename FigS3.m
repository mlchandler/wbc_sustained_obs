% Mitchell Chandler, SIO
% Last updated: 10/03/2022

load ix21_variability
load px30_variability
load px40_variability

%% Running mean
%3-month filter of time series (used for computing annual cycle)
w_size = 3; %window size
ix21_wbc_transport_filt3 = movmean(ix21_wbc_transport_raw,w_size,'EndPoints','fill');
px30_wbc_transport_filt3 = movmean(px30_wbc_transport_raw,w_size,'EndPoints','fill');
px40_wbc_transport_filt3 = movmean(px40_wbc_transport_raw,w_size,'EndPoints','fill');

%% Trends
%Linear trends; monthly time step (dt=1) and 95% CI (alpha=0.05).
[ix21_wbc_transport_coeff,ix21_wbc_transport_trend,ix21_wbc_transport_CI] = linear_trend(time_monthly,ix21_wbc_transport_raw,1,0.05);
[px30_wbc_transport_coeff,px30_wbc_transport_trend,px30_wbc_transport_CI] = linear_trend(time_monthly,px30_wbc_transport_raw,1,0.05);
[px40_wbc_transport_coeff,px40_wbc_transport_trend,px40_wbc_transport_CI] = linear_trend(time_monthly,px40_wbc_transport_raw,1,0.05);

%trends are not significant if the gradient includes 0 within its CI
trend_sig = array2table(NaN(3,1),'VariableNames',{'Trend Significant'},'RowNames',{'Agulhas','EAC','Kuroshio'});
trend_sig{'Agulhas',1} = abs(ix21_wbc_transport_coeff(1)) - abs(ix21_wbc_transport_CI) > 0;
trend_sig{'EAC',1} = abs(px30_wbc_transport_coeff(1)) - abs(px30_wbc_transport_CI) > 0;
trend_sig{'Kuroshio',1} = abs(px40_wbc_transport_coeff(1)) - abs(px40_wbc_transport_CI) > 0;

trend_sig %so trend is only significant for PX40

%compute trend and CI (per year) for Kuroshio
[px40_wbc_transport_coeff(1)*365.25, px40_wbc_transport_CI*365.25]

%% Compute 95% CI bounds for Kuroshio trend
%Upper bound
upper_coeff = [px40_wbc_transport_coeff(1)+px40_wbc_transport_CI, px40_wbc_transport_coeff(2)];
px40_trend_upper_CI = polyval(upper_coeff,time_monthly);
px40_trend_upper_CI = px40_trend_upper_CI-px40_trend_upper_CI(1)+px40_wbc_transport_trend(1);

%Lower bound
lower_coeff = [px40_wbc_transport_coeff(1)-px40_wbc_transport_CI, px40_wbc_transport_coeff(2)];
px40_trend_lower_CI = polyval(lower_coeff,time_monthly);
px40_trend_lower_CI = px40_trend_lower_CI-px40_trend_lower_CI(1)+px40_wbc_transport_trend(1);

%% -- Plot --
fsize = 15;

figure('color','w')
clf

subplot(3,1,1) %ix21
hold on
%time series
plot(time_monthly,ix21_wbc_transport_raw,'k','LineWidth',3)
plot(time_monthly,ix21_wbc_transport_filt3,'Color',rgb('light blue'),'LineWidth',3)
%labels and axis
datetick('x') 
ylim([-90 -10])
yticks(-90:20:-10)
ylabel('Transport [Sv]')
text(datenum('01-Jan-2020'),-3,'(a) Agulhas Current','HorizontalAlignment','right','FontSize',fsize)
box on
grid on
set(gca,'FontSize',fsize)

subplot(3,1,2) %px30
hold on
%time series
plot(time_monthly,px30_wbc_transport_raw,'k','LineWidth',3)
plot(time_monthly,px30_wbc_transport_filt3,'Color',rgb('light blue'),'LineWidth',3)
%labels and axis
datetick('x')
yticks(-40:10:0)
ylabel('Transport [Sv]')
text(datenum('01-Jan-2020'),3,'(b) East Australian Current Current','HorizontalAlignment','right','FontSize',fsize)
box on
grid on
set(gca,'FontSize',fsize)

subplot(3,1,3) %px40
hold on
%time series
plot(time_monthly,px40_wbc_transport_raw,'k','LineWidth',3)
plot(time_monthly,px40_wbc_transport_filt3,'Color',rgb('light blue'),'LineWidth',3)
%trend
plot(time_monthly,px40_wbc_transport_trend,'r','LineWidth',2)
% plot(time_monthly,px40_trend_upper_CI,'r:','LineWidth',1)
% plot(time_monthly,px40_trend_lower_CI,'r:','LineWidth',1)
fill([time_monthly;flipud(time_monthly)],[px40_trend_lower_CI;flipud(px40_trend_upper_CI)],'r','linestyle','none','facealpha',0.15)
%labels and axis
datetick('x')
ylim([10 100])
yticks(10:15:100)
ylabel('Transport [Sv]')
text(datenum('01-Jan-2020'),107,'(c) Kuroshio','HorizontalAlignment','right','FontSize',fsize)
box on
grid on
set(gca,'FontSize',fsize)

