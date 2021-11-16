% Mitchell Chandler, SIO
% Last updated: 08/02/2021

%Compute cumulative transport across a  transect taking the
%depth-integrated velocity measurement to be at the mid-point of each
%coordinate bin (i.e. using a midpoint riemann sum).

%depth_int_vel has dimensions [long/lat] x [time]

function [cumulated_transport,transport_long,transport_lat] = calc_cumulative_transport(long,lat,depth_int_vel,time)
if range(lat) > range(long) %meridional transect
    half_space = mean(diff(lat))/2;
    transport_lat = [lat-half_space; lat(end)+half_space];
    transport_long = interp1(lat,long,transport_lat,'linear','extrap');
elseif range(long) > range(lat) %zonal transect
    half_space = mean(diff(long))/2;
    transport_long = [long-half_space; long(end)+half_space];
    transport_lat = interp1(long,lat,transport_long,'linear','extrap');
end
dist = gsw_distance(transport_long,transport_lat); %m

cumulated_transport = NaN(length(transport_lat),length(time));

for t=1:length(time)
    %Transport for each distance bin along the transect
    transport = 0*transport_lat; %The first transport is 0
    transport(2:end) = depth_int_vel(:,t).*dist;
    
    %Sum transports for accumulated transport from W-E or S-N
    cumulated_transport(:,t) = cumsum(transport,'omitnan')/1E6; %Sv
end
end
