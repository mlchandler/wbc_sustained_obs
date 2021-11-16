% Mitchell Chandler, SIO
% Last updated: 08/02/2021

%Compute depth-integrated velocity along a transect taking the velocity
%measurement to be at the mid-point of each depth bin (i.e. using a
%midpoint riemann sum).

%vel has dimensions [depth] x [lat/long] x [time]

function [depth_int_vel] = calc_depth_int_vel(coord,depth,vel,time)
depth_int_vel = NaN(length(coord),length(time)); %m^2/s
delta_z = diff(depth); 
for t=1:length(time)
    for i=1:length(coord)
        vel_bins = NaN(length(depth),1);
        
        vel_bins(1) = vel(1,i,t) * (diff([0 depth(1)])+delta_z(1)/2); %top depth bin (nb. first Argo pressure level is at 2.5 dbar)
        for z=2:length(depth)-1
           vel_bins(z) =  vel(z,i,t) * (delta_z(z)/2+delta_z(z-1)/2);
        end
        vel_bins(end) = vel(end,i,t) * delta_z(end)/2; %bottom depth bin
        
        depth_int_vel(i,t) = nansum(vel_bins);
    end 
    %set depth-integrated velocity to NaN where velocity is NaN at all depths
    depth_int_vel(find(all(isnan(vel(:,:,t)),1)),t) = NaN; 
end
end