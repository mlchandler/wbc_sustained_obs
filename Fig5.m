% Mitchell Chandler, SIO
% Last updated: 16/06/2022

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

%% Smooth KEI using boxcar filter
%boxcar filter gives the best visual match to smoothing applied by Qiu et al. 2020

w_size = 53; %53 weeks
W = (1/sum(boxcar(w_size)))*boxcar(w_size); 
[KEI_lowpass] = conv_filt(KEI,W,w_size);

%Interpolate KEI to smoothed XAA times (monthly averaging gives results almost the same)
smooth_time = time_monthly(7:end-6);
KEI_interp = interp1(KEI_date,KEI_lowpass,smooth_time);

%% Average velocity cross-sections over +ve and -ve KEI
%Interpolate KEI to unsmoothed XAA times (monthly averaging gives results almost the same)
KEI_composite_interp = interp1(KEI_date,KEI_lowpass,time_monthly);

% idx_pos = find(KEI_composite_interp > 0);
% idx_neg = find(KEI_composite_interp < 0);
idx_pos = find(KEI_composite_interp > std(KEI_composite_interp)/2);
idx_neg = find(KEI_composite_interp < -std(KEI_composite_interp)/2);

idx_long = find(px40_long_nom <= 148);

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

%initialise arrays
store_diff = NaN(size(gvel_diff,1),size(gvel_diff,2),iterations); %[depth x long x iterations]

pos_transport = NaN(iterations,1);
neg_transport = NaN(iterations,1);

pos_max_vel = NaN(iterations,1);
neg_max_vel = NaN(iterations,1);

for i=1:iterations
    %sample with replacement from +ve KEI state to build subsample of same size
    A = datasample(idx_pos,length(idx_pos));
    mean_A = mean(px40_gvel_LKM(:,idx_long,A),3);
        
    %sample with replacement from -ve KEI state to build subsample of same size
    B = datasample(idx_neg,length(idx_neg));
    mean_B = mean(px40_gvel_LKM(:,idx_long,B),3);
    
    %compute difference in velocity between subsampled cross-sections
    store_diff(:,:,i) = mean_A - mean_B;
    
    %compute mean wbc transport of subsamples
    pos_transport(i) = mean(px40_wbc_transport_raw(A));
    neg_transport(i) = mean(px40_wbc_transport_raw(B));
    
    %find max poleward velocity from subsampled composites
    pos_max_vel(i) = max(mean_A,[],'all');
    neg_max_vel(i) = max(mean_B,[],'all');
end

%Find upper and lower percentile values to give CI:
%velocity cross-section composites
lower_bound = prctile(store_diff,alpha/2,3);
upper_bound = prctile(store_diff,100-alpha/2,3);

%Composite difference is significant if within the 95% CI, and the 95% CI does not cross 0.
CI_mask = double(lower_bound.*upper_bound < 0 | gvel_diff < lower_bound | gvel_diff > upper_bound); %(not significant = 1)
% if the CI lower bound and upper bound have different signs then multiplying them will give a negative number
CI_mask(isnan(gvel_diff)) = NaN; %mask NaNs

%wbc transport composites
pos_transport_upper = prctile(pos_transport,100-alpha/2);
pos_transport_lower = prctile(pos_transport,alpha/2);
neg_transport_upper = prctile(neg_transport,100-alpha/2);
neg_transport_lower = prctile(neg_transport,alpha/2);

[pos_transport_upper, neg_transport_upper;
    mean(px40_wbc_transport_raw(idx_pos)), mean(px40_wbc_transport_raw(idx_neg));
    pos_transport_lower, neg_transport_lower]

%maximum poleward composite velocity
pos_vel_upper = prctile(pos_max_vel,100-alpha/2);
pos_vel_lower = prctile(pos_max_vel,alpha/2);
neg_vel_upper = prctile(neg_max_vel,100-alpha/2);
neg_vel_lower = prctile(neg_max_vel,alpha/2);

[pos_vel_upper, neg_vel_upper;
    max(gvel_kei_pos,[],'all'), max(gvel_kei_neg,[],'all');
    pos_vel_lower, neg_vel_lower]

%% -- Plot --
fsize = 15;
xticklength = 6E-3;
cticklength = 0.05;

figure('color','w')
clf

mv = 1;
%+ve KEI composite
subplot(3,1,1) 
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
text(px40_long_nom(idx_long(end)),-1975,'\textbf{(a)} KEI $> \frac{1}{2}\sigma$','interpreter','latex',...
    'FontSize',fsize,'HorizontalAlignment','right','VerticalAlignment','bottom','BackgroundColor','w')
XX = [140:1:148];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-1900 -1300 -700 -100];
yticks(YY)
YT = compose('%.0f m',-YY);
yticklabels(YT)
ylabel('Depth [m]')
set(gca,'TickDir','out','FontSize',fsize)
ax=gca;
ax.XAxis.TickLength = [xticklength xticklength];

%-ve KEI composite
subplot(3,1,2)  
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
text(px40_long_nom(idx_long(end)),-1975,'\textbf{(b)} KEI $< -\frac{1}{2}\sigma$','interpreter','latex',...
    'FontSize',fsize,'HorizontalAlignment','right','VerticalAlignment','bottom','BackgroundColor','w')
XX = [140:1:148];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-1900 -1300 -700 -100];
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
subplot(3,1,3) 
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
c.TickLabels = {'-0.5','','','','','0','','','','','0.5'}; %mv=0.5
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
text(px40_long_nom(idx_long(end)),-1975,'\textbf{(c)} Diff.','interpreter','latex',...
    'FontSize',fsize,'HorizontalAlignment','right','VerticalAlignment','bottom','BackgroundColor','w')
XX = [140:1:148];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
YY = [-1900 -1300 -700 -100];
yticks(YY)
YT = compose('%.0f m',-YY);
yticklabels(YT)
ylabel('Depth [m]')
set(gca,'TickDir','out','FontSize',fsize)
ax=gca;
ax.XAxis.TickLength = [xticklength xticklength];




