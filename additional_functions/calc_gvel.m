% Mitchell Chandler, SIO
% Last updated: 07/02/2021

%Calculate geostrophic velocity based on the thermal wind relation relative
%to a LNM. Derivatives are calculated using centred differencing, with the
%end DH values repeated as ghost points. Velocities below the bathymetry
%are masked.

%v(p) = 1/f * d/dx[H'(p) - H'(p0)]
%u(p) = -1/f * d/dy[H'(p) - H'(p0)]

function [gvel] = calc_gvel(DH,p,long,lat,LNM)
%% Calculate geostrophic velocity
%initialise variables
gvel = NaN*DH;
%distance between points along transect
dist = gsw_distance(long,lat);
%calculate geostrophic velocity using centred differencing for derivatives
for i=1:length(p) %pressure levels
    %outermost DH values are used here as being repeated as ghost points at the end points
    f = gsw_f(lat(1));
    DH_grad = (DH(i,2)-DH(LNM(1),2)) - (DH(i,1)-DH(LNM(1),1)); %use difference in DH between the pressure level and the reference pressure at start location 
    dd = 2*dist(1); %multiply distance by 2 to represent DH being repeated as a ghost point
    gvel(i,1) = 1/f*DH_grad/dd;
    for j=2:length(long)-1
        f = gsw_f(lat(j));
        DH_grad = (DH(i,j+1)-DH(LNM(j),j+1)) - (DH(i,j-1)-DH(LNM(j),j-1)); %use difference in DH between the pressure level and the reference pressure at location j 
        dd = dist(j-1)+dist(j);
        gvel(i,j) = 1/f*DH_grad/dd;
    end
    f = gsw_f(lat(end));
    DH_grad = (DH(i,end)-DH(LNM(end),end)) - (DH(i,end-1)-DH(LNM(end),end-1)); %use difference in DH between the pressure level and the reference pressure at end location
    dd = 2*dist(end); %multiply distance by 2 to represent DH being repeated as a ghost point
    gvel(i,end) = 1/f*DH_grad/dd;
end
%multiply by -1 for meridional transects
if range(lat) > range(long) %meridional 
    gvel = -gvel;
end
%mask velocity values below bathymetry with NaNs
for i=1:length(long) 
    gvel(LNM(i)+1:end,i) = NaN; 
    %the +1 index means that velocities go to the bottom, and that the
    %bathymetry mask better matches a plot of the topography
end
end

