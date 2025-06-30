% WELL LOGS ANALYSIS OF 31/5-7 EOS WELL FOR CO2
% see https://www.equinor.com/news/archive/20201019-sharing-data-northern-lights
% Coded by Ufuk Durmus on 12/2024

close all; clear; clc;

% File Directory
file_path = '.\31_5-7-EOS';

%% Import well logs data
% importing the well logs
measured_logs = importdata([file_path,'\13.Petrophysical_Data_Evaluations\CPI\measured_logs_edited.txt']);
data = measured_logs.data;

% importing processed petrophysical logs
petrophysics_logs = importdata([file_path,'\13.Petrophysical_Data_Evaluations\CPI\petrophysics_logs_edited.txt']);
data2 = petrophysics_logs.data;

% importing core logs
core_logs = importdata([file_path,'\11.Core_Data\core_logs_edited.txt']);
data3 = core_logs;

%% DEFINE VARIABLES
% limit the data such that all logs have the same lenght with petrophysical logs depth
ind1 = 12418; ind2 = 18869;

depth = data(ind1:ind2,1);
cali = data(ind1:ind2,2); cali = clip(cali,0,max(cali));
drho = data(ind1:ind2,3); drho = clip(drho,-1,max(drho));
gr = data(ind1:ind2,6); gr = clip(gr,0,max(gr));
nphi = data(ind1:ind2,8); nphi = clip(nphi,-0.15,max(nphi));
vp = data(ind1:ind2,4); vp = clip(vp,0,max(vp));
rhob = data(ind1:ind2,11); rhob = clip(rhob,0,max(rhob));
vs = data(ind1:ind2,5); vs = clip(vs,0,max(vs));
phit = data2(:,5); phit = clip(phit,0,max(phit));
v_clay = data2(:,7); v_clay = clip(v_clay,0,max(v_clay)); % defined as vshale
% sw = data(:,0); % no need and not found, it is an aquifer
depth2 = data2(:,1); % depth range of petrophysical logs are not the same

% conversion (if necessary)
vp = ((1./vp).*10^6).*0.3048; % m/s
vs = ((1./vs).*10^6).*0.3048; % m/s

v_quartz = 1-v_clay;
VpVs = vp./vs;

% tops
intra_drake = 2585; % md in meters
lower_drake1 = 2610; % from Thompson etal, 2024 ARMA paper
cook4 = 2638; % md in meters, very thin
cook2 = 2642; % md in meters
burton = 2695; % from Thompson and Griffiths, 2024 ARMA paper
johansen4 = 2702; % md in meters
johansen2 = 2752; % md in meters
amundsen = 2818; % md in meters

%% Core measurements from Eos well as a log
depth_core_well = data3(:,1); % in meters
rho_grain_core_well = data3(:,2); % in g/cc
phit_core_well = data3(:,3); phit_core_well = phit_core_well./100; % in v/v
rho_w_core = 1.057; % from PVT report
rhob_core_well = rho_grain_core_well.*(1-phit_core_well) + phit_core_well.*rho_w_core;

%% Display the well logs with home-made function
% classical well logs view for reference
ylimit = [2580 inf];
well_logs_disp_final(depth,depth2,depth_core_well,ylimit,cali,gr,vp,vs,rhob,rhob_core_well,drho,phit,phit_core_well,v_clay,intra_drake,lower_drake1,cook4,cook2,burton,johansen4,johansen2,amundsen)
%% Save the data
% save('well_data_outputs\well_rp_outputs.mat')
