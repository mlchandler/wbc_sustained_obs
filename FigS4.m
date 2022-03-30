% Mitchell Chandler, SIO
% Last updated: 30/03/2022

load px40_velocity
gvel_mean = nanmean(px40_gvel_LKM,3);

%% Find velocity trend
vel_LKM_trend = NaN(length(argo_depth),length(px40_long_nom));
vel_LKM_trend_sig = NaN(length(argo_depth),length(px40_long_nom));
for x=1:length(px40_long_nom)
    for z=1:length(argo_depth)
        Y = squeeze(px40_gvel_LKM(z,x,:));
        [p,Yhat,CI] = linear_trend(time_monthly,Y,1,0.05);
        vel_LKM_trend(z,x) = p(1);
        %only save trend if it is significant
        if abs(p(1)) - abs(CI) > 0 %trend is significant
            vel_LKM_trend_sig(z,x) = p(1);
        elseif isnan(p(1)) %retain NaN values (bathymetry)
            vel_LKM_trend_sig(z,x) = NaN;
            vel_LKM_trend(z,x) = NaN;
        else %trend is not significant
            vel_LKM_trend_sig(z,x) = 0;
        end
    end
end

%convert trends from per day to per year
vel_LKM_trend_sig = vel_LKM_trend_sig*365.25;
vel_LKM_trend = vel_LKM_trend*365.25;

%% Plot
fsize = 11;

xx = [px40_long_nom(1) 143];

C1 = brewermap(256,'*PuOr');
C1(123:134,:) = repmat([1 1 1],12,1);

C2 = brewermap(256,'*RdBu');
C2(123:134,:) = repmat([1 1 1],12,1);


figure('color','w')
clf

%Mean velocity
subplot(2,1,1)
imagesc(px40_long_nom,argo_depth,gvel_mean)
%colours:
colormap(gca,C1)
caxis([-1 1])
c1 = colorbar;
c1.Label.String = 'Velocity [m/s]';
c1.FontSize = fsize;
c1.Label.FontSize = fsize;
%x-axis:
xlim(xx)
XX = [140:1:143];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
%y-axis:
YY = [100:600:1900];
yticks(YY)
YT = compose('%.0f m',YY);
yticklabels(YT)
ylabel('Depth [m]')
%formatting:
set(gca,'FontSize',fsize)
%text:
text(142.75,1800,'(a)','FontSize',fsize,'FontWeight','bold')

%Velocity trend
subplot(2,1,2)
imagesc(px40_long_nom,argo_depth,vel_LKM_trend_sig)
%colours:
colormap(gca,C2)
caxis([-0.1 0.1])
c2 = colorbar;
c2.Label.String = 'Velocity trend [m/s/yr]';
c2.FontSize = fsize;
c2.Label.FontSize = fsize;
%x-axis:
xlim(xx)
XX = [140:1:143];
xticks(XX)
XT = compose('%.0f\\circE',XX);
xticklabels(XT)
%y-axis:
YY = [100:600:1900];
yticks(YY)
YT = compose('%.0f m',YY);
yticklabels(YT)
ylabel('Depth [m]')
%formatting:
set(gca,'FontSize',fsize)
%text:
text(142.75,1800,'(b)','FontSize',fsize,'FontWeight','bold')


