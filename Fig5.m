% Mitchell Chandler, SIO
% Last updated: 09/03/2022

%Colours from Paul Tol (https://personal.sron.nl/~pault/) and Brewermap to ensure colourblind friendly palettes 

load px40_fig2.mat
load px40_variability.mat
load px40_velocity.mat

%% Read in KEI
A = readtable('KE_index_2020_07.csv');
KEI_time = A{:,'Date'};
KEI_date = datenum(KEI_time);
KEI = A{:,'KEI'};
KEIsmoothed = A{:,'KEIsmoothed'};

%% Smooth KEI and Core longitude using boxcar filter
%boxcar filter gives the best visual match to smoothing applied by Qiu et al. 2020

w_size = 53; %53 weeks
W = (1/sum(boxcar(w_size)))*boxcar(w_size); 
[KEI_lowpass] = conv_filt(KEI,W,w_size);

w_size = 13; %13 months
W = (1/sum(boxcar(w_size)))*boxcar(w_size);
[px40_core_long_lowpass] = conv_filt(px40_core_long,W,w_size);
px40_core_long_lowpass = px40_core_long_lowpass(7:end-6);
smooth_time = time_monthly(7:end-6);

%Interpolate KEI to smoothed XAA times (monthly averaging gives results almost the same)
KEI_interp = interp1(KEI_date,KEI_lowpass,smooth_time);

%% Correlate KEI and PX40 core longitude
[r,p,edof] = corr_pval(KEI_interp,px40_core_long_lowpass);
px40_rpd = [r;p;edof] 

%% Correlate KEI and other PX40 WBC metrics
% w_size = 13; %13 months
% W = (1/sum(boxcar(w_size)))*boxcar(w_size);
% % px40_metric = px40_wbc_transport_raw;
% % px40_metric = px40_core_speed_raw;
% % px40_metric = px40_offshore_dev_raw;
% [px40_metric_lowpass] = conv_filt(px40_metric,W,w_size);
% px40_metric_lowpass = px40_metric_lowpass(7:end-6);
% [r,p] = corr_pval(KEI_interp,px40_metric_lowpass)

%% Average velocity cross-sections over +ve and -ve KEI
%Interpolate KEI to unsmoothed XAA times (monthly averaging gives results almost the same)
KEI_composite_interp = interp1(KEI_date,KEI_lowpass,time_monthly);

% idx_pos = find(KEI_composite_interp > 0);
% idx_neg = find(KEI_composite_interp < 0);
idx_pos = find(KEI_composite_interp > std(KEI_composite_interp)/2);
idx_neg = find(KEI_composite_interp < -std(KEI_composite_interp)/2);

idx_long = find(px40_long_nom < 147.5);

gvel_kei_pos = mean(px40_gvel_LKM(:,idx_long,idx_pos),3);
gvel_kei_neg = mean(px40_gvel_LKM(:,idx_long,idx_neg),3);
gvel_diff = gvel_kei_pos - gvel_kei_neg;

%% Use bootstrapping to test whether the mean composites are significantly different
%This method treats each variability state as a separate population 
%and bootstraps from each population to produce samples A and B then
%computes the difference between them. The mean and CI for the bootstrapped
%differences is computed. The difference between the actual composites is
%determined to be significant if it falls within the bootstrapped CI, and
%if the boostrapped CI does not cross 0.

iterations = 1000; %1E5 | 1000
alpha = 0.05; %e.g. 0.1 = 90% CI, 0.05 = 95% CI

store_diff = NaN(size(gvel_diff,1),size(gvel_diff,2),iterations); %[depth x long x iterations]

for i=1:iterations
    %sample with replacement from +ve KEI state to build subsample of same size
    A = datasample(idx_pos,length(idx_pos));
    mean_A = mean(px40_gvel_LKM(:,idx_long,A),3);
    
    %sample with replacement from -ve KEI state to build subsample of same size
    B = datasample(idx_neg,length(idx_neg));
    mean_B = mean(px40_gvel_LKM(:,idx_long,B),3);
    
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
fsize = 17;
xticklength = 4E-3;
cticklength = 0.05;

figure('color','w')
clf

%KEI
subaxis(4,1,1)
yyaxis left
hold on
plot(smooth_time,px40_core_long_lowpass-mean(px40_core_long_lowpass),'Color',[0 153 136]/256,'LineWidth',3)
set(gca,'YColor',[0 153 136]/256)
ylim([-0.65 0.65])
yline(0,'k','LineWidth',1.5)
ylabel({'Deviations in','core longitude [\circE]'})
yyaxis right
hold on
plot(smooth_time,KEI_interp,'Color',[238 119 51]/256,'LineWidth',3)
set(gca,'YColor',[238 119 51]/256)
yline(std(KEI_composite_interp)/2,'--','LineWidth',1.5) %std liines
yline(-std(KEI_composite_interp)/2,'--','LineWidth',1.5)
ylim([-2.2 2.2])
ylabel('KEI')
yline(0,'k')
xticks(datenum(2004,1:12:204,1))
xlim([datenum('01-Jan-2004') datenum('01-Jan-2020')])
datetick('x','keepticks','keeplimits')
box on
text(datenum('01-Sep-2017'),-1.7,'\textbf{(a)}','FontSize',fsize,'interpreter','latex')
set(gca,'XGrid','on','FontSize',fsize)

mv = 1;
%+ve KEI composite
subaxis(4,1,2) 
hold on
set(gca,'Color','k') %background colour (and therefore NaNs)
contourf(px40_long_nom(idx_long),-argo_depth,gvel_kei_pos,-mv:0.1:mv,'LineColor','none') 
contour(px40_long_nom(idx_long),-argo_depth,gvel_kei_pos,[0 0],'LineColor',[0 0 0],'LineWidth',2) %0 m/s
contour(px40_long_nom(idx_long),-argo_depth,gvel_kei_pos,-mv:0.1:-0.1,'--','LineColor',[0 0 0]+0.6,'LineWidth',2) %-ve m/s
contour(px40_long_nom(idx_long),-argo_depth,gvel_kei_pos,0.1:0.1:mv,'LineColor',[0 0 0]+0.6,'LineWidth',2) %+ve m/s
c = colorbar;
caxis([-mv mv])
c.Ticks = [-mv:0.1:mv];
c.TickLabels = {'-1','','','','','','','','','','0','','','','','','','','','','1'}; %mv=1
c.TickLength = cticklength;
ylabel(c,'Velocity [m/s]')
box off
%labels and axis
text(146.6,-1700,'\textbf{(b)} KEI $> \frac{1}{2}\sigma$','FontSize',fsize,'interpreter','latex')
XX = [140:1:147];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-1500 -1000 -500];
yticks(YY)
YT = compose('%.0f m',-YY);
yticklabels(YT)
ylabel('Depth [m]')
set(gca,'TickDir','out','FontSize',fsize)
ax=gca;
ax.XAxis.TickLength = [xticklength xticklength];

%-ve KEI composite
subaxis(4,1,3)  
hold on
set(gca,'Color','k') %background colour (and therefore NaNs)
contourf(px40_long_nom(idx_long),-argo_depth,gvel_kei_neg,-mv:0.1:mv,'LineColor','none')
contour(px40_long_nom(idx_long),-argo_depth,gvel_kei_neg,[0 0],'LineColor',[0 0 0],'LineWidth',2) %0 m/s
contour(px40_long_nom(idx_long),-argo_depth,gvel_kei_neg,-mv:0.1:-0.1,'--','LineColor',[0 0 0]+0.6,'LineWidth',2) %-ve m/s
contour(px40_long_nom(idx_long),-argo_depth,gvel_kei_neg,0.1:0.1:mv,'LineColor',[0 0 0]+0.6,'LineWidth',2) %+ve m/s
c = colorbar;
caxis([-mv mv])
c.Ticks = [-mv:0.1:mv];
c.TickLabels = {'-1','','','','','','','','','','0','','','','','','','','','','1'}; %mv=1
c.TickLength = cticklength;
ylabel(c,'Velocity [m/s]')
box off
%labels and axis
text(146.6,-1700,'\textbf{(c)} KEI $< -\frac{1}{2}\sigma$','FontSize',fsize,'interpreter','latex')
XX = [140:1:147];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-1500 -1000 -500];
yticks(YY)
YT = compose('%.0f m',-YY);
yticklabels(YT)
ylabel('Depth [m]')
set(gca,'TickDir','out','FontSize',fsize)
ax=gca;
ax.XAxis.TickLength = [xticklength xticklength];
%colourmap for above 2 figures
colormap(brewermap(20,'*RdBu'))

%composite difference
subaxis(4,1,4) 
mv = ceil(max(abs(gvel_diff),[],'all')*10)/10;
hold on
set(gca,'Color','k') %background colour (and therefore NaNs)
contourf(px40_long_nom(idx_long),-argo_depth,gvel_diff,256,'LineColor','none')
contour(px40_long_nom(idx_long),-argo_depth,gvel_diff,[0 0],'LineColor',[0 0 0],'LineWidth',2) %0 m/s
contour(px40_long_nom(idx_long),-argo_depth,gvel_diff,-mv:0.1:-0.1,'--','LineColor',[0 0 0]+0.6,'LineWidth',2) %-ve m/s
contour(px40_long_nom(idx_long),-argo_depth,gvel_diff,0.1:0.1:mv,'LineColor',[0 0 0]+0.6,'LineWidth',2) %+ve m/s
colormap(gca,brewermap(256,'*BrBG'))
c = colorbar;
caxis([-mv mv])
c.Ticks = [-mv:0.1:mv];
c.TickLabels = {'-0.6','','','','','','0','','','','','','0.6'}; %mv=0.6
c.TickLength = cticklength;
ylabel(c,'Velocity [m/s]')
box off
%hatching indicates non-significance
[c2,h2] = contour(px40_long_nom(idx_long),-argo_depth,CI_mask,[1 1]);
set(h2,'linestyle','none','Tag','HatchingRegion');
ax1 = gca;
hp = findobj(ax1,'Tag','HatchingRegion');
hh = hatchfill2(hp,'single','HatchColor','k','HatchDensity',20);
%labels and axis
text(146.6,-1700,'\textbf{(d)} Diff.','FontSize',fsize,'interpreter','latex')
XX = [140:1:147];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-1500 -1000 -500];
yticks(YY)
YT = compose('%.0f m',-YY);
yticklabels(YT)
ylabel('Depth [m]')
set(gca,'TickDir','out','FontSize',fsize)
ax=gca;
ax.XAxis.TickLength = [xticklength xticklength];




