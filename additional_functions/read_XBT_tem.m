% Mitchell Chandler, SIO
% Last updated: 26/01/2021

function [XBT_long_recent,XBT_lat_recent,XBT_depth,XBT_time_recent,XBT_tem_recent]...
    = read_XBT_tem(file)
%read in XBT data
XBT_long = ncread(file,'LONGITUDE');
XBT_lat = ncread(file,'LATITUDE');
XBT_depth = ncread(file,'DEPTH');
XBT_time_initial = ncread(file,'TIME');
XBT_tem = ncread(file,'TEM');
%time is days since 01 Jan 1986
XBT_time = XBT_time_initial + datenum('01-Jan-1986');
%Only want data from 2004 onwards
start_2004 = find(XBT_time >= datenum('01-Jan-2004'));
XBT_time_recent = XBT_time(start_2004);
XBT_tem_recent = XBT_tem(:,:,start_2004);
test = size(XBT_lat); %test if constant latitude or longitude file
if test(1)==1 || test(2)==1 %constant latitude (meridional transect)
    XBT_long_recent = XBT_long(:,start_2004);
    XBT_lat_recent = XBT_lat;
else %constant longitude (zonal transect)
    XBT_lat_recent = XBT_lat(:,start_2004);
    XBT_long_recent = XBT_long;
end
end

