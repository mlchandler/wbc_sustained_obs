% Mitchell Chandler, SIO
% Last updated: 25/03/2022

%% Read in data
A = readtable('KE_index_2020_07.csv');
KEI_time = A{:,'Date'};
KEI_date = datenum(KEI_time);
KEI = A{:,'KEI'};

load px40_fig2

%% Smooth KEI using boxcar filter
%boxcar filter gives the best visual match to smoothing applied by Qiu et al. 2020
w_size = 53; %53 weeks
W = (1/sum(boxcar(w_size)))*boxcar(w_size); 
[KEI_lowpass] = conv_filt(KEI,W,w_size);

%% Interpolate KEI to XAA times 
%(monthly averaging gives results almost the same)
KEI_composite_interp = interp1(KEI_date,KEI_lowpass,time_monthly);

%% Plot
fsize = 13;

right_text = datetime('01-Aug-2019');

figure('color','w')
clf
hold on
plot(KEI_time,KEI_lowpass,'k','LineWidth',3)
yline(0,'k','LineWidth',1)
yline(std(KEI_composite_interp)/2,'--','LineWidth',2)
yline(-std(KEI_composite_interp)/2,'--','LineWidth',2)
ylim([-2 2])
ylabel('Kuroshio Extension Index')
xlim([datetime('01-Jan-2004') datetime('01-Jan-2020')])
box on
set(gca,'FontSize',fsize)
text(right_text,1.9,'stable dyanmic state','FontSize',fsize,'Color',[0 0 0]+0.5,...
    'HorizontalAlignment','Right','VerticalAlignment','top')
text(right_text,-1.9,'unstable dyanmic state','FontSize',fsize,'Color',[0 0 0]+0.5,...
    'HorizontalAlignment','Right','VerticalAlignment','bottom')


