% Mitchell Chandler, SIO
% Last updated: 16/09/2021

%Read in and process the Argo trajectories.

%Use the variables for the core Argo floats which use the transmitted gps or
%Argos positions when there are less than 6 surface fixes and the
%extrapolated Argos positions when there are 6+ surface fixes:
%ctsubjday_i_core
%ctsubjday_f_core
%ctsublon_i_core
%ctsublon_f_core
%ctsublat_i_core
%ctsublat_f_core
%csubv_core
%csubvx_core
%csubvy_core
%ctsubpres_core


function [traj_time,traj_long,traj_lat,traj_p,traj_speed,traj_uvel,traj_vvel] = process_argo_traj_v2()
%% Read in Argo trajectories
load argo_traj_Nathalie_16sept2021.mat

%% Take mid-points for lat, long
traj_long = (ctsublon_i_core + ctsublon_f_core)/2;
traj_lat = (ctsublat_i_core + ctsublat_f_core)/2;

%% Take mid-point for time
traj_time = (ctsubjday_i_core + ctsubjday_f_core)/2 + datenum('1-Jan-1950'); %time is days since 01-Jan-1950

%Restrict to 2004--2019 period of interest
time_idx = find(traj_time >= datenum('01-Jan-2004') & traj_time < datenum('01-Jan-2020'));

traj_time = traj_time(time_idx);
traj_long = traj_long(time_idx);
traj_lat = traj_lat(time_idx);

ctsubpres_core_2 = ctsubpres_core(time_idx);
csubv_core_2 = csubv_core(time_idx);
csubvx_core_2 = csubvx_core(time_idx);
csubvy_core_2 = csubvy_core(time_idx);

%% Restrict to pressures between 900 and 1100
%Zilberman et al. 2018: Argo float trajectories are used only if the floats park at pressure between 900 and 1100 dbar.
p_idx = find(ctsubpres_core_2 >=900 & ctsubpres_core_2 <=1100);

traj_time = traj_time(p_idx);
traj_long = traj_long(p_idx);
traj_lat = traj_lat(p_idx);
traj_p = ctsubpres_core_2(p_idx);

csubv_core_2 = csubv_core_2(p_idx);
csubvx_core_2 = csubvx_core_2(p_idx);
csubvy_core_2 = csubvy_core_2(p_idx);

%% Convert velocities from cm/s to m/s
csubv_core_2 = csubv_core_2/100;
csubvx_core_2 = csubvx_core_2/100;
csubvy_core_2 = csubvy_core_2/100;

%consider speeds > 2 m/s to be errors and remove
vel_idx = find(csubv_core_2 < 2);

traj_time = traj_time(vel_idx);
traj_long = traj_long(vel_idx);
traj_lat = traj_lat(vel_idx);
traj_p = traj_p(vel_idx);
traj_speed = csubv_core_2(vel_idx);
traj_uvel = csubvx_core_2(vel_idx);
traj_vvel = csubvy_core_2(vel_idx);
end