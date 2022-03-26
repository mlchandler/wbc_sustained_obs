% Mitchell Chandler, SIO
% Last updated: 25/03/2022

%% Read in data
load ix21_traj_dates
load px30_traj_dates
load px40_traj_dates

unq_yrs = 2004:2019;

%% IX21
%Find year and month of all trajectory data points
[yy21,mm21,~] = ymd(datetime(ix21_traj_time,'ConvertFrom','datenum'));

%Number of profiles each year
ix21_interannual = NaN(length(unq_yrs),length(ix21_bin_long));
for i=1:length(ix21_bin_long)
    for t=1:length(unq_yrs)
        ix21_interannual(t,i) = numel(find(yy21(:,i) == unq_yrs(t)));
    end
end

%Number of profiles each month
ix21_month = NaN(12,length(ix21_bin_long));
for i=1:length(ix21_bin_long)
    for t=1:12
        ix21_month(t,i) = numel(find(mm21(:,i) == t));
    end
end

%% PX30
%Find year and month of all trajectory data points
[yy30,mm30,~] = ymd(datetime(px30_traj_time,'ConvertFrom','datenum'));

%Number of profiles each year
px30_interannual = NaN(length(unq_yrs),length(px30_bin_long));
for i=1:length(px30_bin_long)
    for t=1:length(unq_yrs)
        px30_interannual(t,i) = numel(find(yy30(:,i) == unq_yrs(t)));
    end
end

%Number of profiles each month
px30_month = NaN(12,length(px30_bin_long));
for i=1:length(px30_bin_long)
    for t=1:12
        px30_month(t,i) = numel(find(mm30(:,i) == t));
    end
end

%% PX40
%Find year and month of all trajectory data points
[yy40,mm40,~] = ymd(datetime(px40_traj_time,'ConvertFrom','datenum'));

%Number of profiles each year
px40_interannual = NaN(length(unq_yrs),length(px40_bin_long));
for i=1:length(px40_bin_long)
    for t=1:length(unq_yrs)
        px40_interannual(t,i) = numel(find(yy40(:,i) == unq_yrs(t)));
    end
end

%Number of profiles each month
px40_month = NaN(12,length(px40_bin_long));
for i=1:length(px40_bin_long)
    for t=1:12
        px40_month(t,i) = numel(find(mm40(:,i) == t));
    end
end

%% Plot
cmap1 = brewermap(30,'greens');
cmap1(1,:) = [1 1 1];
cmap1(2:5,:) = [0 0 0;0 0 0;0 0 0;0 0 0]+0.9;

cmap2 = brewermap(20,'blues');
cmap2(1,:) = [1 1 1];
cmap2(2:5,:) = [0 0 0;0 0 0;0 0 0;0 0 0]+0.9;

tlength = 1E-2;
fsize = 13;

figure('color','w')
clf

%ix21 interannual
subaxis(2,3,1)
imagesc(ix21_bin_long,unq_yrs,ix21_interannual)
caxis([0 30])
colormap(gca,cmap1)
%y-axis
yticks(unq_yrs)
yticklabels({'2004','','2006','','2008','','2010','','2012','','2014','','2016','','2018',''})
ylabel('Year')
%x-axis
XX = [30:5:55];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
%manually draw box to avoid top and right ticks
box off
yline(2019.5)
xline(ix21_bin_long(end)+0.25)
%format
set(gca,'YDir','normal','TickDir','out','TickLength',[tlength tlength],'FontSize',fsize)
%text
text(ix21_bin_long(end)-0.08,2019.3,'(a)','FontSize',fsize,...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'BackgroundColor','k','Color','w','FontWeight','bold')
title('IX21','FontWeight','normal')

%px30 interannual
subaxis(2,3,2)
imagesc(px30_bin_long,unq_yrs,px30_interannual)
caxis([0 30])
colormap(gca,cmap1)
%y-axis
yticks(unq_yrs)
yticklabels({'2004','','2006','','2008','','2010','','2012','','2014','','2016','','2018',''})
%x-axis
XX = [155:5:175];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
%manually draw box to avoid top and right ticks
box off
yline(2019.5)
xline(px30_bin_long(end)+0.25)
%format
set(gca,'YDir','normal','TickDir','out','TickLength',[tlength tlength],'FontSize',fsize)
%text
text(px30_bin_long(end)-0.08,2019.3,'(b)','FontSize',fsize,...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'BackgroundColor','k','Color','w','FontWeight','bold')
title('PX30','FontWeight','normal')

%px40 interannual
subaxis(2,3,3)
imagesc(px40_bin_long,unq_yrs,px40_interannual)
caxis([0 30])
colormap(gca,cmap1)
%y-axis
yticks(unq_yrs)
yticklabels({'2004','','2006','','2008','','2010','','2012','','2014','','2016','','2018',''})
%x-axis
XX = [140:15:200];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
%manually draw box to avoid top and right ticks
box off
yline(2019.5)
xline(px40_bin_long(end)+0.25)
%format
set(gca,'YDir','normal','TickDir','out','TickLength',[tlength tlength],'FontSize',fsize)
%text
text(px40_bin_long(end)-0.5,2019.3,'(c)','FontSize',fsize,...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'BackgroundColor','k','Color','w','FontWeight','bold')
title('PX40','FontWeight','normal')
%colourbar
c1 = colorbar;
c1.Label.String = 'Number of Argo trajectories';
c1.TickLength = tlength;
c1.FontSize = fsize;
c1.Label.FontSize = fsize;
c1.TickDirection = 'out';
set(c1,'Position', [0.91, 0.525, 0.015, 0.375]);


%ix21 annual
subaxis(2,3,4)
imagesc(ix21_bin_long,1:12,ix21_month)
caxis([0 20])
colormap(gca,cmap2)
%y-axis
yticks(1:12)
yticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel('Month')
%x-axis
XX = [30:5:55];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
%manually draw box to avoid top and right ticks
box off
yline(12.5)
xline(ix21_bin_long(end)+0.25)
%format
set(gca,'YDir','normal','TickDir','out','TickLength',[tlength tlength],'FontSize',fsize)
%text
text(ix21_bin_long(end)-0.08,12.3,'(d)','FontSize',fsize,...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'BackgroundColor','k','Color','w','FontWeight','bold')

%px30 annual
subaxis(2,3,5)
imagesc(px30_bin_long,1:12,px30_month)
caxis([0 20])
colormap(gca,cmap2)
%y-axis
yticks(1:12)
yticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
%x-axis
XX = [155:5:175];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
%manually draw box to avoid top and right ticks
box off
yline(12.5)
xline(px30_bin_long(end)+0.25)
%format
set(gca,'YDir','normal','TickDir','out','TickLength',[tlength tlength],'FontSize',fsize)
%text
text(px30_bin_long(end)-0.08,12.3,'(e)','FontSize',fsize,...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'BackgroundColor','k','Color','w','FontWeight','bold')

%px40 annual
subaxis(2,3,6)
imagesc(px40_bin_long,1:12,px40_month)
caxis([0 20])
colormap(gca,cmap2)
%y-axis
yticks(1:12)
yticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
%x-axis
XX = [140:15:200];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
%manually draw box to avoid top and right ticks
box off
yline(12.5)
xline(px40_bin_long(end)+0.25)
%format
set(gca,'YDir','normal','TickDir','out','TickLength',[tlength tlength],'FontSize',fsize)
%text
text(px40_bin_long(end)-0.5,12.3,'(f)','FontSize',fsize,...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'BackgroundColor','k','Color','w','FontWeight','bold')
%colourbar
c2 = colorbar;
c2.Label.String = 'Number of Argo trajectories';
c2.TickLength = tlength;
c2.FontSize = fsize;
c2.Label.FontSize = fsize;
c2.TickDirection = 'out';
set(c2,'Position', [0.91, 0.1, 0.015, 0.375]);


