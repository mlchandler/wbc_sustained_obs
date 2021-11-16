% Mitchell Chandler, SIO
% Last updated: 02/08/2021

%Estimate the (sampling) error on the mean WBC transport using the standard
%error of the Argo sub-surface velocity in each bin out to the mean
%offshore edge of the WBC.

function [mean_transport,transport_uncertainty,transport_STD,transport_SE] = ...
    mean_error_transport(wbc_transport,argo_sub_SE,argo_sub_long,long_nom,lat_nom,inshore_long,offshore_long)

%Interpolate argo sub-surface velocity SE to transect longitude points using nearest neighbour interpolation
argo_sub_SE_nearest = interp1(argo_sub_long,argo_sub_SE,long_nom,'nearest','extrap');

%Restrict to between inshore edge and mean WBC offshore edge
argo_sub_SE_wbc = argo_sub_SE_nearest(long_nom >= inshore_long & long_nom <= offshore_long);
long_nom_wbc = long_nom(long_nom >= inshore_long & long_nom <= offshore_long);
lat_nom_wbc = lat_nom(long_nom >= inshore_long & long_nom <= offshore_long);

%Read in bathymetry
topo_long = ncread('etopo6.nc','LON6');
topo_lat = ncread('etopo6.nc','LAT6101_1700');
topo_bath = ncread('etopo6.nc','BATH6');

%Depth at each point in WBC
topo_wbc = -1*interp2(topo_lat,topo_long,topo_bath,lat_nom_wbc,long_nom_wbc);
topo_wbc(topo_wbc > 1975) = 1975; %maximum depth of 1975-m
    
figure()
yyaxis left
plot(long_nom_wbc,-topo_wbc,'b','LineWidth',2)
ylabel('Bathymetry depth [m]')
yyaxis right
hold on 
plot(argo_sub_long,argo_sub_SE,'-o')
plot(long_nom_wbc,argo_sub_SE_wbc,'r-','LineWidth',2)
ylabel('Argo subsurface velocity SE  [m/s]')
xlabel('Longitude [\circE]')
xlim([long_nom(1) long_nom(50)])
xline(offshore_long,'--')
xline(inshore_long,'--')
box on

%Velocity SE constant in depth so multiply by depth to depth-integrate
argo_sub_SE_wbcDI = argo_sub_SE_wbc.*topo_wbc;

%Transport uncertainty
[cumulated_transport,transport_long,transport_lat] = calc_cumulative_transport(long_nom_wbc,lat_nom_wbc,argo_sub_SE_wbcDI,0);

%Outputs:
mean_transport = mean(wbc_transport); %Mean transport
transport_uncertainty = cumulated_transport(end); %Transport uncertainty
transport_STD = std(wbc_transport); %Transport standard deviation

%Find SE of mean transport for comparison with assessed uncertainty
[~,~,eDOF] = corr_pval(wbc_transport,wbc_transport); %take into account auto-correlation
transport_SE = std(wbc_transport)/sqrt(eDOF); %Transport SE
end
