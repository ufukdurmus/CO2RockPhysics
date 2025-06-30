% ROCK MECHANICS DATA ANALYSIS OF 31/5-7 EOS WELL FOR CO2
% download the data here https://co2datashare.org/dataset/rock-mechanics-data-for-31-5-7-eos-well-core
% Coded by Ufuk Durmus on 8/2024

close all; clear; clc;

% File Directory
file_path = 'C:\yourpath\data';

%% Import rock mechanics data
% Brazilian test for shale and sandstones
core_brazil = importdata([file_path,'\Brazil_tests.xlsx']);
core_brazil_data = core_brazil.data;

% Shale samples
shale_caust_T2528 = importdata([file_path,'\Shale_CAUST\31_5-7_T2528_CAUST-data.xlsx']);
shale_caust_T2528_data = shale_caust_T2528.data.TimeSeries;

shale_ciu_T2492 = importdata([file_path,'\Shale_CIU\31_5-7_T2492_CIU-data.xlsx']);
shale_ciu_T2492_data = shale_ciu_T2492.data.TimeSeries;

shale_ciu_T2514 = importdata([file_path,'\Shale_CIU\31_5-7_T2514_CIU-data.xlsx']);
shale_ciu_T2514_data = shale_ciu_T2514.data.TimeSeries;

shale_ciu_T2522 = importdata([file_path,'\Shale_CIU\31_5-7_T2522_CIU-data.xlsx']);
shale_ciu_T2522_data = shale_ciu_T2522.data.TimeSeries;

% NOTE: T2529 and T2533 samples are horizontal
shale_ciu_T2529 = importdata([file_path,'\Shale_CIU\31_5-7_T2529_CIU-data.xlsx']);
shale_ciu_T2529_data = shale_ciu_T2529.data.TimeSeries;

shale_ciu_T2533 = importdata([file_path,'\Shale_CIU\31_5-7_T2533_CIU-data.xlsx']);
shale_ciu_T2533_data = shale_ciu_T2533.data.TimeSeries;

% Sandstone samples
ss_caust_T2505 = importdata([file_path,'\SS_CAUST\31_5-7_T2505_UST-data_wCalc.xlsx']);
ss_caust_T2505_data = ss_caust_T2505.data.TimeSeries;

ss_caust_T2507 = importdata([file_path,'\SS_CAUST\31_5-7_T2507_UST-data_wCalc.xlsx']);
ss_caust_T2507_data = ss_caust_T2507.data.TimeSeries;

ss_caust_T2511 = importdata([file_path,'\SS_CAUST\31_5-7_T2511_UST-data_wCalc.xlsx']);
ss_caust_T2511_data = ss_caust_T2511.data.TimeSeries;

ss_cid_T2461 = importdata([file_path,'\SS_CID\31_5-7_T2461_CID-data.xlsx']);
ss_cid_T2461_data = ss_cid_T2461.data.TimeSeries;

ss_cid_T2463 = importdata([file_path,'\SS_CID\31_5-7_T2463_CID-data.xlsx']);
ss_cid_T2463_data = ss_cid_T2463.data.TimeSeries;

ss_cid_T2464 = importdata([file_path,'\SS_CID\31_5-7_T2464_CID-data.xlsx']);
ss_cid_T2464_data = ss_cid_T2464.data.TimeSeries;

ss_cid_T2466 = importdata([file_path,'\SS_CID\31_5-7_T2466_CID-data.xlsx']);
ss_cid_T2466_data = ss_cid_T2466.data.TimeSeries;

ss_cid_T2467 = importdata([file_path,'\SS_CID\31_5-7_T2467_CID-data.xlsx']);
ss_cid_T2467_data = ss_cid_T2467.data.TimeSeries;

ss_cid_T2469 = importdata([file_path,'\SS_CID\31_5-7_T2469_CID-data.xlsx']);
ss_cid_T2469_data = ss_cid_T2469.data.TimeSeries;

ss_cid_T2470 = importdata([file_path,'\SS_CID\31_5-7_T2470_CID-data.xlsx']);
ss_cid_T2470_data = ss_cid_T2470.data.TimeSeries;

ss_cid_T2472 = importdata([file_path,'\SS_CID\31_5-7_T2472_CID-data.xlsx']);
ss_cid_T2472_data = ss_cid_T2472.data.TimeSeries;

ss_cid_T2476 = importdata([file_path,'\SS_CID\31_5-7_T2476_CID-data.xlsx']);
ss_cid_T2476_data = ss_cid_T2476.data.TimeSeries;

ss_hydro_T2479 = importdata([file_path,'\SS_Hydrostat\31_5-7_T2479_Hydrostatic-data.xlsx']);
ss_hydro_T2479_data = ss_hydro_T2479.data.TimeSeries;

ss_hydro_T2480 = importdata([file_path,'\SS_Hydrostat\31_5-7_T2480_Hydrostatic-data.xlsx']);
ss_hydro_T2480_data = ss_hydro_T2480.data.TimeSeries;

ss_hydro_T2482 = importdata([file_path,'\SS_Hydrostat\31_5-7_T2482_Hydrostatic-data.xlsx']);
ss_hydro_T2482_data = ss_hydro_T2482.data.TimeSeries;

%% CORE DENSITY ESTIMATION
% unit weight of solids in kN/m3
unit_weight_solid = [shale_caust_T2528.data.SingleInputValues(6,6)...
    shale_ciu_T2492.data.SingleInputValues(6,6)...
    shale_ciu_T2514.data.SingleInputValues(6,6)...
    shale_ciu_T2522.data.SingleInputValues(6,6)...
    shale_ciu_T2529.data.SingleInputValues(6,6)...
    shale_ciu_T2533.data.SingleInputValues(6,6)...
    ss_caust_T2505.data.SingleInputValues(6,6)...
    ss_caust_T2507.data.SingleInputValues(6,6)...
    ss_caust_T2511.data.SingleInputValues(6,6)...
    ss_cid_T2461.data.SingleInputValues(6,6)...
    ss_cid_T2463.data.SingleInputValues(6,6)...
    ss_cid_T2464.data.SingleInputValues(6,6)...
    ss_cid_T2466.data.SingleInputValues(6,6)...
    ss_cid_T2467.data.SingleInputValues(6,6)...
    ss_cid_T2469.data.SingleInputValues(6,6)...
    ss_cid_T2470.data.SingleInputValues(6,6)...
    ss_cid_T2472.data.SingleInputValues(6,6)...
    ss_cid_T2476.data.SingleInputValues(6,6)...
    ss_hydro_T2479.data.SingleInputValues(6,6)...
    ss_hydro_T2480.data.SingleInputValues(6,6)...
    ss_hydro_T2482.data.SingleInputValues(6,6)];

% unit weight of particles in kN/m3
unit_weight_particles = [shale_caust_T2528.data.SingleInputValues(10,1)...
    shale_ciu_T2492.data.SingleInputValues(10,1)...
    shale_ciu_T2514.data.SingleInputValues(10,1)...
    shale_ciu_T2522.data.SingleInputValues(10,1)...
    shale_ciu_T2529.data.SingleInputValues(10,1)...
    shale_ciu_T2533.data.SingleInputValues(10,1)...
    ss_caust_T2505.data.SingleInputValues(10,1)...
    ss_caust_T2507.data.SingleInputValues(10,1)...
    ss_caust_T2511.data.SingleInputValues(10,1)...
    ss_cid_T2461.data.SingleInputValues(10,1)...
    ss_cid_T2463.data.SingleInputValues(10,1)...
    ss_cid_T2464.data.SingleInputValues(10,1)...
    ss_cid_T2466.data.SingleInputValues(10,1)...
    ss_cid_T2467.data.SingleInputValues(10,1)...
    ss_cid_T2469.data.SingleInputValues(10,1)...
    ss_cid_T2470.data.SingleInputValues(10,1)...
    ss_cid_T2472.data.SingleInputValues(10,1)...
    ss_cid_T2476.data.SingleInputValues(10,1)...
    ss_hydro_T2479.data.SingleInputValues(10,1)...
    ss_hydro_T2480.data.SingleInputValues(10,1)...
    ss_hydro_T2482.data.SingleInputValues(10,1)];

% Convert unit weight to density from kN/m3 to g/cc
rho_solid =((unit_weight_solid.*1000)./9.80665)./1000;
rho_particles =((unit_weight_particles.*1000)./9.80665)./1000;

% initial porosity of the samples in percentage
por_initial = [shale_caust_T2528.data.SingleInputValues(7,6)...
    shale_ciu_T2492.data.SingleInputValues(7,6)...
    shale_ciu_T2514.data.SingleInputValues(7,6)...
    shale_ciu_T2522.data.SingleInputValues(7,6)...
    shale_ciu_T2529.data.SingleInputValues(7,6)...
    shale_ciu_T2533.data.SingleInputValues(7,6)...
    ss_caust_T2505.data.SingleInputValues(7,6)...
    ss_caust_T2507.data.SingleInputValues(7,6)...
    ss_caust_T2511.data.SingleInputValues(7,6)...
    ss_cid_T2461.data.SingleInputValues(7,6)...
    ss_cid_T2463.data.SingleInputValues(7,6)...
    ss_cid_T2464.data.SingleInputValues(7,6)...
    ss_cid_T2466.data.SingleInputValues(7,6)...
    ss_cid_T2467.data.SingleInputValues(7,6)...
    ss_cid_T2469.data.SingleInputValues(7,6)...
    ss_cid_T2470.data.SingleInputValues(7,6)...
    ss_cid_T2472.data.SingleInputValues(7,6)...
    ss_cid_T2476.data.SingleInputValues(7,6)...
    ss_hydro_T2479.data.SingleInputValues(7,6)...
    ss_hydro_T2480.data.SingleInputValues(7,6)...
    ss_hydro_T2482.data.SingleInputValues(7,6)];

phi_initial = por_initial./100;

% Sw of the samples (sandstones are drained and shales are fully saturated)
sw_core = zeros(1,21); sw_core(1,1:6) = 1; sw_core(1,7:21) = 0;

% density of water (for shales 120K ppm, for SS 80K ppm salinity data from excels)
rho_water = zeros(1,21); rho_water(1,1:6) = 1.092; rho_water(1,7:21) = 1.06;

% estimate bulk density of the samples (both shales and ss 21 total)
rho_dry = rho_particles.*(1-phi_initial);
rhob = rho_dry + phi_initial.*sw_core.*rho_water; % assumption from the report shales are saturated, sands are dry

%% STRESS-STRAIN
% SHALES
% ev = vertical strain
ev_shale_caust_T2528 = shale_caust_T2528_data(:,4);
% eh = horizontal strain
eh_shale_caust_T2528 = shale_caust_T2528_data(:,5);
% evol = volumetric strain
evol_shale_caust_T2528 = shale_caust_T2528_data(:,6);
% sv = vertical stress
sv_shale_caust_T2528 = shale_caust_T2528_data(:,9);
% sh = horizontal stress
sh_shale_caust_T2528 = shale_caust_T2528_data(:,10);
% sshear = shear stress
sshear_shale_caust_T2528 = shale_caust_T2528_data(:,11);
% % differential stress
% sdiff_shale_caust_T2528 = shale_caust_T2528_data(:,9)-shale_caust_T2528_data(:,10);

% find the data that makes sense (with figures) and limit the data number
ev_shale_ciu_T2492 = shale_ciu_T2492_data((12091:end),4);
eh_shale_ciu_T2492 = shale_ciu_T2492_data((12091:end),5);
evol_shale_ciu_T2492 = shale_ciu_T2492_data((12091:end),6);
sv_shale_ciu_T2492 = shale_ciu_T2492_data((12091:end),9);
sh_shale_ciu_T2492 = shale_ciu_T2492_data((12091:end),10);
sshear_shale_ciu_T2492 = shale_ciu_T2492_data((12091:end),11);

ev_shale_ciu_T2514 = shale_ciu_T2514_data((8770:end),4);
eh_shale_ciu_T2514 = shale_ciu_T2514_data((8770:end),5);
evol_shale_ciu_T2514 = shale_ciu_T2514_data((8770:end),6);
sv_shale_ciu_T2514 = shale_ciu_T2514_data((8770:end),9);
sh_shale_ciu_T2514 = shale_ciu_T2514_data((8770:end),10);
sshear_shale_ciu_T2514 = shale_ciu_T2514_data((8770:end),11);

ev_shale_ciu_T2522 = shale_ciu_T2522_data((4000:end),4);
eh_shale_ciu_T2522 = shale_ciu_T2522_data((4000:end),5);
evol_shale_ciu_T2522 = shale_ciu_T2522_data((4000:end),6);
sv_shale_ciu_T2522 = shale_ciu_T2522_data((4000:end),9);
sh_shale_ciu_T2522 = shale_ciu_T2522_data((4000:end),10);
sshear_shale_ciu_T2522 = shale_ciu_T2522_data((4000:end),11);

ev_shale_ciu_T2529 = shale_ciu_T2529_data((4500:end),4);
eh_shale_ciu_T2529 = shale_ciu_T2529_data((4500:end),5);
evol_shale_ciu_T2529 = shale_ciu_T2529_data((4500:end),6);
sv_shale_ciu_T2529 = shale_ciu_T2529_data((4500:end),9);
sh_shale_ciu_T2529 = shale_ciu_T2529_data((4500:end),10);
sshear_shale_ciu_T2529 = shale_ciu_T2529_data((4500:end),11);

ev_shale_ciu_T2533 = shale_ciu_T2533_data((3400:end),4);
eh_shale_ciu_T2533 = shale_ciu_T2533_data((3400:end),5); % weird data, see the curve
evol_shale_ciu_T2533 = shale_ciu_T2533_data((3400:end),6);
sv_shale_ciu_T2533 = shale_ciu_T2533_data((3400:end),9);
sh_shale_ciu_T2533 = shale_ciu_T2533_data((3400:end),10);
sshear_shale_ciu_T2533 = shale_ciu_T2533_data((3400:end),11);

% SANDSTONES
ev_ss_caust_T2505 = ss_caust_T2505_data(:,4);
eh_ss_caust_T2505 = ss_caust_T2505_data(:,5);
evol_ss_caust_T2505 = ss_caust_T2505_data(:,6);
sv_ss_caust_T2505 = ss_caust_T2505_data(:,9);
sh_ss_caust_T2505 = ss_caust_T2505_data(:,10);
sshear_ss_caust_T2505 = ss_caust_T2505_data(:,11);

ev_ss_caust_T2507 = ss_caust_T2507_data(:,4);
eh_ss_caust_T2507 = ss_caust_T2507_data(:,5);
evol_ss_caust_T2507 = ss_caust_T2507_data(:,6);
sv_ss_caust_T2507 = ss_caust_T2507_data(:,9);
sh_ss_caust_T2507 = ss_caust_T2507_data(:,10);
sshear_ss_caust_T2507 = ss_caust_T2507_data(:,11);

ev_ss_caust_T2511 = ss_caust_T2511_data(:,4);
eh_ss_caust_T2511 = ss_caust_T2511_data(:,5);
evol_ss_caust_T2511 = ss_caust_T2511_data(:,6);
sv_ss_caust_T2511 = ss_caust_T2511_data(:,9);
sh_ss_caust_T2511 = ss_caust_T2511_data(:,10);
sshear_ss_caust_T2511 = ss_caust_T2511_data(:,11);

ev_ss_cid_T2461 = ss_cid_T2461_data((2400:end),4);
eh_ss_cid_T2461 = ss_cid_T2461_data((2400:end),5);
evol_ss_cid_T2461 = ss_cid_T2461_data((2400:end),6);
sv_ss_cid_T2461 = ss_cid_T2461_data((2400:end),9);
sh_ss_cid_T2461 = ss_cid_T2461_data((2400:end),10);
sshear_ss_cid_T2461 = ss_cid_T2461_data((2400:end),11);

ev_ss_cid_T2463 = ss_cid_T2463_data((900:end),4);
eh_ss_cid_T2463 = ss_cid_T2463_data((900:end),5);
evol_ss_cid_T2463 = ss_cid_T2463_data((900:end),6);
sv_ss_cid_T2463 = ss_cid_T2463_data((900:end),9);
sh_ss_cid_T2463 = ss_cid_T2463_data((900:end),10);
sshear_ss_cid_T2463 = ss_cid_T2463_data((900:end),11);

ev_ss_cid_T2464 = ss_cid_T2464_data((6320:end),4);
eh_ss_cid_T2464 = ss_cid_T2464_data((6320:end),5);
evol_ss_cid_T2464 = ss_cid_T2464_data((6320:end),6);
sv_ss_cid_T2464 = ss_cid_T2464_data((6320:end),9);
sh_ss_cid_T2464 = ss_cid_T2464_data((6320:end),10);
sshear_ss_cid_T2464 = ss_cid_T2464_data((6320:end),11);

ev_ss_cid_T2466 = ss_cid_T2466_data((1331:end),4);
eh_ss_cid_T2466 = ss_cid_T2466_data((1331:end),5);
evol_ss_cid_T2466 = ss_cid_T2466_data((1331:end),6);
sv_ss_cid_T2466 = ss_cid_T2466_data((1331:end),9);
sh_ss_cid_T2466 = ss_cid_T2466_data((1331:end),10);
sshear_ss_cid_T2466 = ss_cid_T2466_data((1331:end),11);

ev_ss_cid_T2467 = ss_cid_T2467_data((651:end),4);
eh_ss_cid_T2467 = ss_cid_T2467_data((651:end),5);
evol_ss_cid_T2467 = ss_cid_T2467_data((651:end),6);
sv_ss_cid_T2467 = ss_cid_T2467_data((651:end),9);
sh_ss_cid_T2467 = ss_cid_T2467_data((651:end),10);
sshear_ss_cid_T2467 = ss_cid_T2467_data((651:end),11);

ev_ss_cid_T2469 = ss_cid_T2469_data((2335:end),4);
eh_ss_cid_T2469 = ss_cid_T2469_data((2335:end),5);
evol_ss_cid_T2469 = ss_cid_T2469_data((2335:end),6);
sv_ss_cid_T2469 = ss_cid_T2469_data((2335:end),9);
sh_ss_cid_T2469 = ss_cid_T2469_data((2335:end),10);
sshear_ss_cid_T2469 = ss_cid_T2469_data((2335:end),11);

ev_ss_cid_T2470 = ss_cid_T2470_data((1219:end),4);
eh_ss_cid_T2470 = ss_cid_T2470_data((1219:end),5);
evol_ss_cid_T2470 = ss_cid_T2470_data((1219:end),6);
sv_ss_cid_T2470 = ss_cid_T2470_data((1219:end),9);
sh_ss_cid_T2470 = ss_cid_T2470_data((1219:end),10);
sshear_ss_cid_T2470 = ss_cid_T2470_data((1219:end),11);

ev_ss_cid_T2472 = ss_cid_T2472_data((1100:end),4);
eh_ss_cid_T2472 = ss_cid_T2472_data((1100:end),5);
evol_ss_cid_T2472 = ss_cid_T2472_data((1100:end),6);
sv_ss_cid_T2472 = ss_cid_T2472_data((1100:end),9);
sh_ss_cid_T2472 = ss_cid_T2472_data((1100:end),10);
sshear_ss_cid_T2472 = ss_cid_T2472_data((1100:end),11);

ev_ss_cid_T2476 = ss_cid_T2476_data((1489:end),4);
eh_ss_cid_T2476 = ss_cid_T2476_data((1489:end),5);
evol_ss_cid_T2476 = ss_cid_T2476_data((1489:end),6);
sv_ss_cid_T2476 = ss_cid_T2476_data((1489:end),9);
sh_ss_cid_T2476 = ss_cid_T2476_data((1489:end),10);
sshear_ss_cid_T2476 = ss_cid_T2476_data((1489:end),11);

ev_ss_hydro_T2479 = ss_hydro_T2479_data(:,4);
eh_ss_hydro_T2479 = ss_hydro_T2479_data(:,5);
evol_ss_hydro_T2479 = ss_hydro_T2479_data(:,6);
sv_ss_hydro_T2479 = ss_hydro_T2479_data(:,9);
sh_ss_hydro_T2479 = ss_hydro_T2479_data(:,10);
sshear_ss_hydro_T2479 = ss_hydro_T2479_data(:,11);

ev_ss_hydro_T2480 = ss_hydro_T2480_data(:,4);
eh_ss_hydro_T2480 = ss_hydro_T2480_data(:,5);
evol_ss_hydro_T2480 = ss_hydro_T2480_data(:,6);
sv_ss_hydro_T2480 = ss_hydro_T2480_data(:,9);
sh_ss_hydro_T2480 = ss_hydro_T2480_data(:,10);
sshear_ss_hydro_T2480 = ss_hydro_T2480_data(:,11);

ev_ss_hydro_T2482 = ss_hydro_T2482_data(:,4);
eh_ss_hydro_T2482 = ss_hydro_T2482_data(:,5);
evol_ss_hydro_T2482 = ss_hydro_T2482_data(:,6);
sv_ss_hydro_T2482 = ss_hydro_T2482_data(:,9);
sh_ss_hydro_T2482 = ss_hydro_T2482_data(:,10);
sshear_ss_hydro_T2482 = ss_hydro_T2482_data(:,11);

%% Core Cij and Anisotropy Estimation from Velocities
% velocities in m/s
% NOTE: limit the data to remove the plastic region measurements, see the raw data curves
% NOTE: velocity plots - clipped the data (see the raw plots)
[rowind_T2528,~] = find(isnan(shale_caust_T2528_data(:,20))); %find the nan indices and then remove from the data
vp_axial_shale_caust_T2528 = shale_caust_T2528_data(:,20); vp_axial_shale_caust_T2528(rowind_T2528,:)=[];
vs_axial_shale_caust_T2528 = shale_caust_T2528_data(:,21); vs_axial_shale_caust_T2528(rowind_T2528,:)=[];
vp_radial_shale_caust_T2528 = shale_caust_T2528_data(:,22); vp_radial_shale_caust_T2528(rowind_T2528,:)=[];
vs_radial_shale_caust_T2528 = shale_caust_T2528_data(:,23); vs_radial_shale_caust_T2528(rowind_T2528,:)=[];
sdiff_shale_caust_T2528 = shale_caust_T2528_data(:,9)-shale_caust_T2528_data(:,10); sdiff_shale_caust_T2528(rowind_T2528,:)=[];

[rowind_T2492,~] = find(isnan(shale_ciu_T2492_data(:,20)));
vp_axial_shale_ciu_T2492 = shale_ciu_T2492_data(:,20); vp_axial_shale_ciu_T2492(rowind_T2492,:)=[]; vp_axial_shale_ciu_T2492=vp_axial_shale_ciu_T2492(2000:2480);
vs_axial_shale_ciu_T2492 = shale_ciu_T2492_data(:,21); vs_axial_shale_ciu_T2492(rowind_T2492,:)=[]; vs_axial_shale_ciu_T2492=vs_axial_shale_ciu_T2492(2000:2480);
sdiff_shale_ciu_T2492 = shale_ciu_T2492_data(:,9)-shale_ciu_T2492_data(:,10); sdiff_shale_ciu_T2492(rowind_T2492,:)=[]; sdiff_shale_ciu_T2492=sdiff_shale_ciu_T2492(2000:2480);

[rowind_T2514,~] = find(isnan(shale_ciu_T2514_data(:,20)));
vp_axial_shale_ciu_T2514 = shale_ciu_T2514_data(:,20); vp_axial_shale_ciu_T2514(rowind_T2514,:)=[]; vp_axial_shale_ciu_T2514=vp_axial_shale_ciu_T2514(800:890);
sdiff_shale_ciu_T2514 = shale_ciu_T2514_data(:,9)-shale_ciu_T2514_data(:,10); sdiff_shale_ciu_T2514(rowind_T2514,:)=[]; sdiff_shale_ciu_T2514=sdiff_shale_ciu_T2514(800:890);

[rowind_T2522,~] = find(isnan(shale_ciu_T2522_data(:,20)));
vp_axial_shale_ciu_T2522 = shale_ciu_T2522_data(:,20); vp_axial_shale_ciu_T2522(rowind_T2522,:)=[]; vp_axial_shale_ciu_T2522=vp_axial_shale_ciu_T2522(179:296);
vs_axial_shale_ciu_T2522 = shale_ciu_T2522_data(:,21); vs_axial_shale_ciu_T2522(rowind_T2522,:)=[]; vs_axial_shale_ciu_T2522=vs_axial_shale_ciu_T2522(179:296);
sdiff_shale_ciu_T2522 = shale_ciu_T2522_data(:,9)-shale_ciu_T2522_data(:,10); sdiff_shale_ciu_T2522(rowind_T2522,:)=[]; sdiff_shale_ciu_T2522=sdiff_shale_ciu_T2522(179:296);

[rowind_T2529,~] = find(isnan(shale_ciu_T2529_data(:,20)));
vp_axial_shale_ciu_T2529 = shale_ciu_T2529_data(:,20); vp_axial_shale_ciu_T2529(rowind_T2529,:)=[]; vp_axial_shale_ciu_T2529=vp_axial_shale_ciu_T2529(474:564);
sdiff_shale_ciu_T2529 = shale_ciu_T2529_data(:,9)-shale_ciu_T2529_data(:,10); sdiff_shale_ciu_T2529(rowind_T2529,:)=[]; sdiff_shale_ciu_T2529=sdiff_shale_ciu_T2529(474:564);

[rowind_T2533,~] = find(isnan(shale_ciu_T2533_data(:,20)));
vp_axial_shale_ciu_T2533 = shale_ciu_T2533_data(:,20); vp_axial_shale_ciu_T2533(rowind_T2533,:)=[]; vp_axial_shale_ciu_T2533=vp_axial_shale_ciu_T2533(472:512);
sdiff_shale_ciu_T2533 = shale_ciu_T2533_data(:,9)-shale_ciu_T2533_data(:,10); sdiff_shale_ciu_T2533(rowind_T2533,:)=[]; sdiff_shale_ciu_T2533=sdiff_shale_ciu_T2533(472:512);

[rowind_T2505,~] = find(isnan(ss_caust_T2505_data(:,20)));
vp_axial_ss_caust_T2505 = ss_caust_T2505_data(:,20); vp_axial_ss_caust_T2505(rowind_T2505,:)=[]; vp_axial_ss_caust_T2505=vp_axial_ss_caust_T2505(640:1348);
vs_axial_ss_caust_T2505 = ss_caust_T2505_data(:,21); vs_axial_ss_caust_T2505(rowind_T2505,:)=[]; vs_axial_ss_caust_T2505=vs_axial_ss_caust_T2505(640:1348);
vp_radial_ss_caust_T2505 = ss_caust_T2505_data(:,22); vp_radial_ss_caust_T2505(rowind_T2505,:)=[]; vp_radial_ss_caust_T2505=vp_radial_ss_caust_T2505(640:1348);
vs_radial_ss_caust_T2505 = ss_caust_T2505_data(:,23); vs_radial_ss_caust_T2505(rowind_T2505,:)=[]; vs_radial_ss_caust_T2505=vs_radial_ss_caust_T2505(640:1348);
sdiff_ss_caust_T2505 = ss_caust_T2505_data(:,9)-ss_caust_T2505_data(:,10); sdiff_ss_caust_T2505(rowind_T2505,:)=[]; sdiff_ss_caust_T2505=sdiff_ss_caust_T2505(640:1348);

[rowind_T2507,~] = find(isnan(ss_caust_T2507_data(:,20)));
vs_axial_ss_caust_T2507 = ss_caust_T2507_data(:,20); vs_axial_ss_caust_T2507(rowind_T2507,:)=[]; vs_axial_ss_caust_T2507=vs_axial_ss_caust_T2507(480:767);
vp_radial_ss_caust_T2507 = ss_caust_T2507_data(:,21); vp_radial_ss_caust_T2507(rowind_T2507,:)=[]; vp_radial_ss_caust_T2507=vp_radial_ss_caust_T2507(480:767);
vs_radial_ss_caust_T2507 = ss_caust_T2507_data(:,22); vs_radial_ss_caust_T2507(rowind_T2507,:)=[]; vs_radial_ss_caust_T2507=vs_radial_ss_caust_T2507(480:767);
sdiff_ss_caust_T2507 = ss_caust_T2507_data(:,9)-ss_caust_T2507_data(:,10); sdiff_ss_caust_T2507(rowind_T2507,:)=[]; sdiff_ss_caust_T2507=sdiff_ss_caust_T2507(480:767);

[rowind_T2511,~] = find(isnan(ss_caust_T2511_data(:,20)));
vp_axial_ss_caust_T2511 = ss_caust_T2511_data(:,20); vp_axial_ss_caust_T2511(rowind_T2511,:)=[]; vp_axial_ss_caust_T2511=vp_axial_ss_caust_T2511(970:1335);
vs_axial_ss_caust_T2511 = ss_caust_T2511_data(:,21); vs_axial_ss_caust_T2511(rowind_T2511,:)=[]; vs_axial_ss_caust_T2511=vs_axial_ss_caust_T2511(970:1335); 
vp_radial_ss_caust_T2511 = ss_caust_T2511_data(:,22); vp_radial_ss_caust_T2511(rowind_T2511,:)=[]; vp_radial_ss_caust_T2511=vp_radial_ss_caust_T2511(970:1335);
vs_radial_ss_caust_T2511 = ss_caust_T2511_data(:,23); vs_radial_ss_caust_T2511(rowind_T2511,:)=[]; vs_radial_ss_caust_T2511=vs_radial_ss_caust_T2511(970:1335);
sdiff_ss_caust_T2511 = ss_caust_T2511_data(:,9)-ss_caust_T2511_data(:,10); sdiff_ss_caust_T2511(rowind_T2511,:)=[]; sdiff_ss_caust_T2511=sdiff_ss_caust_T2511(970:1335);

[rowind_T2461,~] = find(isnan(ss_cid_T2461_data(:,25))); 
vp_axial_ss_cid_T2461 = ss_cid_T2461_data(:,25); vp_axial_ss_cid_T2461(rowind_T2461,:)=[]; vp_axial_ss_cid_T2461=vp_axial_ss_cid_T2461(30:69);
vs_axial_ss_cid_T2461 = ss_cid_T2461_data(:,26); vs_axial_ss_cid_T2461(rowind_T2461,:)=[]; vs_axial_ss_cid_T2461=vs_axial_ss_cid_T2461(30:69);
sdiff_ss_cid_T2461 = ss_cid_T2461_data(:,9)-ss_cid_T2461_data(:,10); sdiff_ss_cid_T2461(rowind_T2461,:)=[]; sdiff_ss_cid_T2461=sdiff_ss_cid_T2461(30:69);

[rowind_T2463,~] = find(isnan(ss_cid_T2463_data(:,20))); 
vp_axial_ss_cid_T2463 = ss_cid_T2463_data(:,20); vp_axial_ss_cid_T2463(rowind_T2463,:)=[]; vp_axial_ss_cid_T2463=vp_axial_ss_cid_T2463(110:130);
vs_axial_ss_cid_T2463 = ss_cid_T2463_data(:,21); vs_axial_ss_cid_T2463(rowind_T2463,:)=[]; vs_axial_ss_cid_T2463=vs_axial_ss_cid_T2463(110:130);
sdiff_ss_cid_T2463 = ss_cid_T2463_data(:,9)-ss_cid_T2463_data(:,10); sdiff_ss_cid_T2463(rowind_T2463,:)=[]; sdiff_ss_cid_T2463=sdiff_ss_cid_T2463(110:130);

[rowind_T2464,~] = find(isnan(ss_cid_T2464_data(:,20))); 
vp_axial_ss_cid_T2464 = ss_cid_T2464_data(:,20); vp_axial_ss_cid_T2464(rowind_T2464,:)=[]; vp_axial_ss_cid_T2464=vp_axial_ss_cid_T2464(198:250);
vs_axial_ss_cid_T2464 = ss_cid_T2464_data(:,21); vs_axial_ss_cid_T2464(rowind_T2464,:)=[]; vs_axial_ss_cid_T2464=vs_axial_ss_cid_T2464(198:250);
sdiff_ss_cid_T2464 = ss_cid_T2464_data(:,9)-ss_cid_T2464_data(:,10); sdiff_ss_cid_T2464(rowind_T2464,:)=[]; sdiff_ss_cid_T2464=sdiff_ss_cid_T2464(198:250);

[rowind_T2466,~] = find(isnan(ss_cid_T2466_data(:,20))); 
vp_axial_ss_cid_T2466 = ss_cid_T2466_data(:,20); vp_axial_ss_cid_T2466(rowind_T2466,:)=[]; vp_axial_ss_cid_T2466=vp_axial_ss_cid_T2466(60:74);
vs_axial_ss_cid_T2466 = ss_cid_T2466_data(:,21); vs_axial_ss_cid_T2466(rowind_T2466,:)=[]; vs_axial_ss_cid_T2466=vs_axial_ss_cid_T2466(60:74);
sdiff_ss_cid_T2466 = ss_cid_T2466_data(:,9)-ss_cid_T2466_data(:,10); sdiff_ss_cid_T2466(rowind_T2466,:)=[]; sdiff_ss_cid_T2466=sdiff_ss_cid_T2466(60:74);

[rowind_T2467,~] = find(isnan(ss_cid_T2467_data(:,20))); 
vp_axial_ss_cid_T2467 = ss_cid_T2467_data(:,20); vp_axial_ss_cid_T2467(rowind_T2467,:)=[]; vp_axial_ss_cid_T2467=vp_axial_ss_cid_T2467(83:104);
vs_axial_ss_cid_T2467 = ss_cid_T2467_data(:,21); vs_axial_ss_cid_T2467(rowind_T2467,:)=[]; vs_axial_ss_cid_T2467=vs_axial_ss_cid_T2467(83:104);
sdiff_ss_cid_T2467 = ss_cid_T2467_data(:,9)-ss_cid_T2467_data(:,10); sdiff_ss_cid_T2467(rowind_T2467,:)=[]; sdiff_ss_cid_T2467=sdiff_ss_cid_T2467(83:104);

[rowind_T2469,~] = find(isnan(ss_cid_T2469_data(:,21))); 
vp_axial_ss_cid_T2469 = ss_cid_T2469_data(:,21); vp_axial_ss_cid_T2469(rowind_T2469,:)=[]; vp_axial_ss_cid_T2469=vp_axial_ss_cid_T2469(135:155);
vs_axial_ss_cid_T2469 = ss_cid_T2469_data(:,22); vs_axial_ss_cid_T2469(rowind_T2469,:)=[]; vs_axial_ss_cid_T2469=vs_axial_ss_cid_T2469(135:155);
sdiff_ss_cid_T2469 = ss_cid_T2469_data(:,9)-ss_cid_T2469_data(:,10); sdiff_ss_cid_T2469(rowind_T2469,:)=[]; sdiff_ss_cid_T2469=sdiff_ss_cid_T2469(135:155);

[rowind_T2470,~] = find(isnan(ss_cid_T2470_data(:,20))); 
vp_axial_ss_cid_T2470 = ss_cid_T2470_data(:,20); vp_axial_ss_cid_T2470(rowind_T2470,:)=[]; vp_axial_ss_cid_T2470=vp_axial_ss_cid_T2470(165:190);
vs_axial_ss_cid_T2470 = ss_cid_T2470_data(:,21); vs_axial_ss_cid_T2470(rowind_T2470,:)=[]; vs_axial_ss_cid_T2470=vs_axial_ss_cid_T2470(165:190);
sdiff_ss_cid_T2470 = ss_cid_T2470_data(:,9)-ss_cid_T2470_data(:,10); sdiff_ss_cid_T2470(rowind_T2470,:)=[]; sdiff_ss_cid_T2470=sdiff_ss_cid_T2470(165:190);

[rowind_T2472,~] = find(isnan(ss_cid_T2472_data(:,20))); 
vp_axial_ss_cid_T2472 = ss_cid_T2472_data(:,20); vp_axial_ss_cid_T2472(rowind_T2472,:)=[]; vp_axial_ss_cid_T2472=vp_axial_ss_cid_T2472(325:375);
vs_axial_ss_cid_T2472 = ss_cid_T2472_data(:,21); vs_axial_ss_cid_T2472(rowind_T2472,:)=[]; vs_axial_ss_cid_T2472=vs_axial_ss_cid_T2472(325:375);
sdiff_ss_cid_T2472 = ss_cid_T2472_data(:,9)-ss_cid_T2472_data(:,10); sdiff_ss_cid_T2472(rowind_T2472,:)=[]; sdiff_ss_cid_T2472=sdiff_ss_cid_T2472(325:375);

[rowind_T2476,~] = find(isnan(ss_cid_T2476_data(:,20))); 
vp_axial_ss_cid_T2476 = ss_cid_T2476_data(:,20); vp_axial_ss_cid_T2476(rowind_T2476,:)=[]; vp_axial_ss_cid_T2476=vp_axial_ss_cid_T2476(40:98);
vs_axial_ss_cid_T2476 = ss_cid_T2476_data(:,21); vs_axial_ss_cid_T2476(rowind_T2476,:)=[]; vs_axial_ss_cid_T2476=vs_axial_ss_cid_T2476(40:98);
sdiff_ss_cid_T2476 = ss_cid_T2476_data(:,9)-ss_cid_T2476_data(:,10); sdiff_ss_cid_T2476(rowind_T2476,:)=[]; sdiff_ss_cid_T2476=sdiff_ss_cid_T2476(40:98);

[rowind_T2479,~] = find(isnan(ss_hydro_T2479_data(:,20))); 
vp_axial_ss_hydro_T2479 = ss_hydro_T2479_data(:,20); vp_axial_ss_hydro_T2479(rowind_T2479,:)=[];
vs_axial_ss_hydro_T2479 = ss_hydro_T2479_data(:,21); vs_axial_ss_hydro_T2479(rowind_T2479,:)=[];
sdiff_ss_hydro_T2479 = ss_hydro_T2479_data(:,9)-0; sdiff_ss_hydro_T2479(rowind_T2479,:)=[];

[rowind_T2480,~] = find(isnan(ss_hydro_T2480_data(:,20))); 
vp_axial_ss_hydro_T2480 = ss_hydro_T2480_data(:,20); vp_axial_ss_hydro_T2480(rowind_T2480,:)=[];
vs_axial_ss_hydro_T2480 = ss_hydro_T2480_data(:,21); vs_axial_ss_hydro_T2480(rowind_T2480,:)=[];
sdiff_ss_hydro_T2480 = ss_hydro_T2480_data(:,9)-0; sdiff_ss_hydro_T2480(rowind_T2480,:)=[];

[rowind_T2482,~] = find(isnan(ss_hydro_T2482_data(:,20))); 
vp_axial_ss_hydro_T2482 = ss_hydro_T2482_data(:,20); vp_axial_ss_hydro_T2482(rowind_T2482,:)=[];
vs_axial_ss_hydro_T2482 = ss_hydro_T2482_data(:,21); vs_axial_ss_hydro_T2482(rowind_T2482,:)=[];
sdiff_ss_hydro_T2482 = ss_hydro_T2482_data(:,9)-0; sdiff_ss_hydro_T2482(rowind_T2482,:)=[];

% Estimate Stiffness coefficients
c33_shale_caust_T2528 = ((vp_axial_shale_caust_T2528./1000).^2).*rhob(1);
c55_shale_caust_T2528 = ((vs_axial_shale_caust_T2528./1000).^2).*rhob(1);
c11_shale_caust_T2528 = ((vp_radial_shale_caust_T2528./1000).^2).*rhob(1);
c66_shale_caust_T2528 = ((vs_radial_shale_caust_T2528./1000).^2).*rhob(1);

c33_shale_ciu_T2492 = ((vp_axial_shale_ciu_T2492./1000).^2).*rhob(2);
c55_shale_ciu_T2492 = ((vs_axial_shale_ciu_T2492./1000).^2).*rhob(2);

c33_shale_ciu_T2514 = ((vp_axial_shale_ciu_T2514./1000).^2).*rhob(3);

c33_shale_ciu_T2522 = ((vp_axial_shale_ciu_T2522./1000).^2).*rhob(4);
c55_shale_ciu_T2522 = ((vs_axial_shale_ciu_T2522./1000).^2).*rhob(4);

c11_shale_ciu_T2529 = ((vp_axial_shale_ciu_T2529./1000).^2).*rhob(5); % hor

c11_shale_ciu_T2533 = ((vp_axial_shale_ciu_T2533./1000).^2).*rhob(6); % hor

c33_ss_caust_T2505 = ((vp_axial_ss_caust_T2505./1000).^2).*rhob(7);
c55_ss_caust_T2505 = ((vs_axial_ss_caust_T2505./1000).^2).*rhob(7);
c11_ss_caust_T2505 = ((vp_radial_ss_caust_T2505./1000).^2).*rhob(7);
c66_ss_caust_T2505 = ((vs_radial_ss_caust_T2505./1000).^2).*rhob(7);

c55_ss_caust_T2507 = ((vs_axial_ss_caust_T2507./1000).^2).*rhob(8);
c11_ss_caust_T2507 = ((vp_radial_ss_caust_T2507./1000).^2).*rhob(8);
c66_ss_caust_T2507 = ((vs_radial_ss_caust_T2507./1000).^2).*rhob(8);

c33_ss_caust_T2511 = ((vp_axial_ss_caust_T2511./1000).^2).*rhob(9);
c55_ss_caust_T2511 = ((vs_axial_ss_caust_T2511./1000).^2).*rhob(9);
c11_ss_caust_T2511 = ((vp_radial_ss_caust_T2511./1000).^2).*rhob(9);
c66_ss_caust_T2511 = ((vs_radial_ss_caust_T2511./1000).^2).*rhob(9);

c33_ss_cid_T2461 = ((vp_axial_ss_cid_T2461./1000).^2).*rhob(10);
c55_ss_cid_T2461 = ((vs_axial_ss_cid_T2461./1000).^2).*rhob(10);

c33_ss_cid_T2463 = ((vp_axial_ss_cid_T2463./1000).^2).*rhob(11);
c55_ss_cid_T2463 = ((vs_axial_ss_cid_T2463./1000).^2).*rhob(11);

c33_ss_cid_T2464 = ((vp_axial_ss_cid_T2464./1000).^2).*rhob(12);
c55_ss_cid_T2464 = ((vs_axial_ss_cid_T2464./1000).^2).*rhob(12);

c33_ss_cid_T2466 = ((vp_axial_ss_cid_T2466./1000).^2).*rhob(13);
c55_ss_cid_T2466 = ((vs_axial_ss_cid_T2466./1000).^2).*rhob(13);

c33_ss_cid_T2467 = ((vp_axial_ss_cid_T2467./1000).^2).*rhob(14);
c55_ss_cid_T2467 = ((vs_axial_ss_cid_T2467./1000).^2).*rhob(14);

c33_ss_cid_T2469 = ((vp_axial_ss_cid_T2469./1000).^2).*rhob(15);
c55_ss_cid_T2469 = ((vs_axial_ss_cid_T2469./1000).^2).*rhob(15);

c33_ss_cid_T2470 = ((vp_axial_ss_cid_T2470./1000).^2).*rhob(16);
c55_ss_cid_T2470 = ((vs_axial_ss_cid_T2470./1000).^2).*rhob(16);

c33_ss_cid_T2472 = ((vp_axial_ss_cid_T2472./1000).^2).*rhob(17);
c55_ss_cid_T2472 = ((vs_axial_ss_cid_T2472./1000).^2).*rhob(17);

c33_ss_cid_T2476 = ((vp_axial_ss_cid_T2476./1000).^2).*rhob(18);
c55_ss_cid_T2476 = ((vs_axial_ss_cid_T2476./1000).^2).*rhob(18);

c33_ss_hydro_T2479 = ((vp_axial_ss_hydro_T2479./1000).^2).*rhob(19);
c55_ss_hydro_T2479 = ((vs_axial_ss_hydro_T2479./1000).^2).*rhob(19);

c33_ss_hydro_T2480 = ((vp_axial_ss_hydro_T2480./1000).^2).*rhob(20);
c55_ss_hydro_T2480 = ((vs_axial_ss_hydro_T2480./1000).^2).*rhob(20);

c33_ss_hydro_T2482 = ((vp_axial_ss_hydro_T2482./1000).^2).*rhob(21);
c55_ss_hydro_T2482 = ((vs_axial_ss_hydro_T2482./1000).^2).*rhob(21);

%% Calculating DYNAMIC Young's Modulus and Poisson's Ratio from Cores

% SANDSTONES & SHALES ISOTROPIC
% Even though they measured axial and radial velocities for SS core samples
% Here I assume SS is mostly isotropic and shales just to see the results
% NOTE that for some samples Vp or Vs is missing and could not compute dynamic E or PR
pr_iso_ss_cid_T2461 = 0.5.*((vp_axial_ss_cid_T2461./vs_axial_ss_cid_T2461).^2-2)./((vp_axial_ss_cid_T2461./vs_axial_ss_cid_T2461).^2-1);
E_iso_ss_cid_T2461 = 2.*c55_ss_cid_T2461.*(1+pr_iso_ss_cid_T2461);

pr_iso_ss_cid_T2463 = 0.5.*((vp_axial_ss_cid_T2463./vs_axial_ss_cid_T2463).^2-2)./((vp_axial_ss_cid_T2463./vs_axial_ss_cid_T2463).^2-1);
E_iso_ss_cid_T2463 = 2.*c55_ss_cid_T2463.*(1+pr_iso_ss_cid_T2463);

pr_iso_ss_cid_T2464 = 0.5.*((vp_axial_ss_cid_T2464./vs_axial_ss_cid_T2464).^2-2)./((vp_axial_ss_cid_T2464./vs_axial_ss_cid_T2464).^2-1);
E_iso_ss_cid_T2464 = 2.*c55_ss_cid_T2464.*(1+pr_iso_ss_cid_T2464);

pr_iso_ss_cid_T2466 = 0.5.*((vp_axial_ss_cid_T2466./vs_axial_ss_cid_T2466).^2-2)./((vp_axial_ss_cid_T2466./vs_axial_ss_cid_T2466).^2-1);
E_iso_ss_cid_T2466 = 2.*c55_ss_cid_T2466.*(1+pr_iso_ss_cid_T2466);

pr_iso_ss_cid_T2467 = 0.5.*((vp_axial_ss_cid_T2467./vs_axial_ss_cid_T2467).^2-2)./((vp_axial_ss_cid_T2467./vs_axial_ss_cid_T2467).^2-1);
E_iso_ss_cid_T2467 = 2.*c55_ss_cid_T2467.*(1+pr_iso_ss_cid_T2467);

pr_iso_ss_cid_T2469 = 0.5.*((vp_axial_ss_cid_T2469./vs_axial_ss_cid_T2469).^2-2)./((vp_axial_ss_cid_T2469./vs_axial_ss_cid_T2469).^2-1);
E_iso_ss_cid_T2469 = 2.*c55_ss_cid_T2469.*(1+pr_iso_ss_cid_T2469);

pr_iso_ss_cid_T2470 = 0.5.*((vp_axial_ss_cid_T2470./vs_axial_ss_cid_T2470).^2-2)./((vp_axial_ss_cid_T2470./vs_axial_ss_cid_T2470).^2-1);
E_iso_ss_cid_T2470 = 2.*c55_ss_cid_T2470.*(1+pr_iso_ss_cid_T2470);

pr_iso_ss_cid_T2472 = 0.5.*((vp_axial_ss_cid_T2472./vs_axial_ss_cid_T2472).^2-2)./((vp_axial_ss_cid_T2472./vs_axial_ss_cid_T2472).^2-1);
E_iso_ss_cid_T2472 = 2.*c55_ss_cid_T2472.*(1+pr_iso_ss_cid_T2472);

pr_iso_ss_cid_T2476 = 0.5.*((vp_axial_ss_cid_T2476./vs_axial_ss_cid_T2476).^2-2)./((vp_axial_ss_cid_T2476./vs_axial_ss_cid_T2476).^2-1);
E_iso_ss_cid_T2476 = 2.*c55_ss_cid_T2476.*(1+pr_iso_ss_cid_T2476);

pr_iso_ss_hydro_T2479 = 0.5.*((vp_axial_ss_hydro_T2479./vs_axial_ss_hydro_T2479).^2-2)./((vp_axial_ss_hydro_T2479./vs_axial_ss_hydro_T2479).^2-1);
E_iso_ss_hydro_T2479 = 2.*c55_ss_hydro_T2479.*(1+pr_iso_ss_hydro_T2479);

pr_iso_ss_hydro_T2480 = 0.5.*((vp_axial_ss_hydro_T2480./vs_axial_ss_hydro_T2480).^2-2)./((vp_axial_ss_hydro_T2480./vs_axial_ss_hydro_T2480).^2-1);
E_iso_ss_hydro_T2480 = 2.*c55_ss_hydro_T2480.*(1+pr_iso_ss_hydro_T2480);

pr_iso_ss_hydro_T2482 = 0.5.*((vp_axial_ss_hydro_T2482./vs_axial_ss_hydro_T2482).^2-2)./((vp_axial_ss_hydro_T2482./vs_axial_ss_hydro_T2482).^2-1);
E_iso_ss_hydro_T2482 = 2.*c55_ss_hydro_T2482.*(1+pr_iso_ss_hydro_T2482);

pr_iso_ss_caust_T2505 = 0.5.*((vp_axial_ss_caust_T2505./vs_axial_ss_caust_T2505).^2-2)./((vp_axial_ss_caust_T2505./vs_axial_ss_caust_T2505).^2-1);
E_iso_ss_caust_T2505 = 2.*c55_ss_caust_T2505.*(1+pr_iso_ss_caust_T2505);

pr_iso_ss_caust_T2511 = 0.5.*((vp_axial_ss_caust_T2511./vs_axial_ss_caust_T2511).^2-2)./((vp_axial_ss_caust_T2511./vs_axial_ss_caust_T2511).^2-1);
E_iso_ss_caust_T2511 = 2.*c55_ss_caust_T2511.*(1+pr_iso_ss_caust_T2511);

pr_iso_shale_ciu_T2492 = 0.5.*((vp_axial_shale_ciu_T2492./vs_axial_shale_ciu_T2492).^2-2)./((vp_axial_shale_ciu_T2492./vs_axial_shale_ciu_T2492).^2-1);
E_iso_shale_ciu_T2492 = 2.*c55_shale_ciu_T2492.*(1+pr_iso_shale_ciu_T2492);

pr_iso_shale_ciu_T2522 = 0.5.*((vp_axial_shale_ciu_T2522./vs_axial_shale_ciu_T2522).^2-2)./((vp_axial_shale_ciu_T2522./vs_axial_shale_ciu_T2522).^2-1);
E_iso_shale_ciu_T2522 = 2.*c55_shale_ciu_T2522.*(1+pr_iso_shale_ciu_T2522);

pr_iso_shale_caust_T2528 = 0.5.*((vp_axial_shale_caust_T2528./vs_axial_shale_caust_T2528).^2-2)./((vp_axial_shale_caust_T2528./vs_axial_shale_caust_T2528).^2-1);
E_iso_shale_caust_T2528 = 2.*c55_shale_caust_T2528.*(1+pr_iso_shale_caust_T2528);

%% ESTIMATION OF STATIC YOUNG'S MODULUS AND POISSON'S RATIO
% SHALES
% find the max value and the indices in the corresponding array
% [peak1,ind_shale_ciu_T2492] = max((sv_shale_ciu_T2492-sh_shale_ciu_T2492));
[peak1,~] = max((sv_shale_ciu_T2492-sh_shale_ciu_T2492));
[peak1,ind1] = min(abs((sv_shale_ciu_T2492-sh_shale_ciu_T2492)-peak1/2));
% define indices for E50 and Etangential = %70-%30 from the peak stress
%ind1 = round(ind1); %i1tan70 = round(ind1*1.4); i1tan30 = round(ind1*0.6);
% estimate two types of young's modulus with differential pressure
E50_shale_ciu_T2492 = (sv_shale_ciu_T2492(ind1)-sh_shale_ciu_T2492(ind1))...
    /(ev_shale_ciu_T2492(ind1));
pr_shale_ciu_T2492 = (-eh_shale_ciu_T2492(ind1))/(ev_shale_ciu_T2492(ind1));
% Etan_shale_ciu_T2492 = ((sv_shale_ciu_T2492(i1tan70)-sh_shale_ciu_T2492(i1tan70))...
%     -(sv_shale_ciu_T2492(i1tan30)-sh_shale_ciu_T2492(i1tan30)))...
%     /(ev_shale_ciu_T2492(i1tan70)-ev_shale_ciu_T2492(i1tan30));

[peak2,~] = max((sv_shale_ciu_T2514-sh_shale_ciu_T2514));
[peak2,ind2] = min(abs((sv_shale_ciu_T2514-sh_shale_ciu_T2514)-peak2/2));
E50_shale_ciu_T2514 = (sv_shale_ciu_T2514(ind2)-sh_shale_ciu_T2514(ind2))...
    /(ev_shale_ciu_T2514(ind2));
pr_shale_ciu_T2514 = (-eh_shale_ciu_T2514(ind2))/(ev_shale_ciu_T2514(ind2));

[peak3,~] = max((sv_shale_ciu_T2522-sh_shale_ciu_T2522));
[peak3,ind3] = min(abs((sv_shale_ciu_T2522-sh_shale_ciu_T2522)-peak3/2));
E50_shale_ciu_T2522 = (sv_shale_ciu_T2522(ind3)-sh_shale_ciu_T2522(ind3))...
    /(ev_shale_ciu_T2522(ind3));
pr_shale_ciu_T2522 = (-eh_shale_ciu_T2522(ind3))/(ev_shale_ciu_T2522(ind3));

[peak4,~] = max((sv_shale_ciu_T2529-sh_shale_ciu_T2529));
[peak4,ind4] = min(abs((sv_shale_ciu_T2529-sh_shale_ciu_T2529)-peak4/2));
ind4 = ind4-800; % normally gives pr > 0.5, changed the indice
E50_shale_ciu_T2529 = (sv_shale_ciu_T2529(ind4)-sh_shale_ciu_T2529(ind4))...
    /(ev_shale_ciu_T2529(ind4));
pr_shale_ciu_T2529 = (-eh_shale_ciu_T2529(ind4))/(ev_shale_ciu_T2529(ind4));

[peak5,~] = max((sv_shale_ciu_T2533-sh_shale_ciu_T2533));
[peak5,ind5] = min(abs((sv_shale_ciu_T2533-sh_shale_ciu_T2533)-peak5/2));
E50_shale_ciu_T2533 = (sv_shale_ciu_T2533(ind5)-sh_shale_ciu_T2533(ind5))...
    /(ev_shale_ciu_T2533(ind5));
pr_shale_ciu_T2533 = (-eh_shale_ciu_T2533(ind5))/(ev_shale_ciu_T2533(ind5));

[peak6,~] = max((sv_ss_cid_T2461-sh_ss_cid_T2461));
[peak6,ind6] = min(abs((sv_ss_cid_T2461-sh_ss_cid_T2461)-peak6/2));
E50_ss_cid_T2461 = (sv_ss_cid_T2461(ind6)-sh_ss_cid_T2461(ind6))...
    /(ev_ss_cid_T2461(ind6));
pr_ss_cid_T2461 = (-eh_ss_cid_T2461(ind6))/(ev_ss_cid_T2461(ind6));

[peak7,~] = max((sv_ss_cid_T2463-sh_ss_cid_T2463));
[peak7,ind7] = min(abs((sv_ss_cid_T2463-sh_ss_cid_T2463)-peak7/2));
E50_ss_cid_T2463 = (sv_ss_cid_T2463(ind7)-sh_ss_cid_T2463(ind7))...
    /(ev_ss_cid_T2463(ind7));
pr_ss_cid_T2463 = (-eh_ss_cid_T2463(ind7))/(ev_ss_cid_T2463(ind7));

[peak8,~] = max((sv_ss_cid_T2464-sh_ss_cid_T2464));
[peak8,ind8] = min(abs((sv_ss_cid_T2464-sh_ss_cid_T2464)-peak8/2));
E50_ss_cid_T2464 = (sv_ss_cid_T2464(ind8)-sh_ss_cid_T2464(ind8))...
    /(ev_ss_cid_T2464(ind8));
pr_ss_cid_T2464 = (-eh_ss_cid_T2464(ind8))/(ev_ss_cid_T2464(ind8));

[peak9,~] = max((sv_ss_cid_T2466-sh_ss_cid_T2466)); 
[peak9,ind9] = min(abs((sv_ss_cid_T2466-sh_ss_cid_T2466)-(peak9/2-0.1))); % added 0.1 see the curve for the reason
E50_ss_cid_T2466 = (sv_ss_cid_T2466(ind9)-sh_ss_cid_T2466(ind9))...
    /(ev_ss_cid_T2466(ind9));
pr_ss_cid_T2466 = (-eh_ss_cid_T2466(ind9))/(ev_ss_cid_T2466(ind9));

[peak10,~] = max((sv_ss_cid_T2467-sh_ss_cid_T2467));
[peak10,ind10] = min(abs((sv_ss_cid_T2467-sh_ss_cid_T2467)-peak10/2));
E50_ss_cid_T2467 = (sv_ss_cid_T2467(ind10)-sh_ss_cid_T2467(ind10))...
    /(ev_ss_cid_T2467(ind10));
pr_ss_cid_T2467 = (-eh_ss_cid_T2467(ind10))/(ev_ss_cid_T2467(ind10));

[peak11,~] = max((sv_ss_cid_T2469-sh_ss_cid_T2469));
[peak11,ind11] = min(abs((sv_ss_cid_T2469-sh_ss_cid_T2469)-peak11/2));
E50_ss_cid_T2469 = (sv_ss_cid_T2469(ind11)-sh_ss_cid_T2469(ind11))...
    /(ev_ss_cid_T2469(ind11));
pr_ss_cid_T2469 = (-eh_ss_cid_T2469(ind11))/(ev_ss_cid_T2469(ind11));

[peak12,~] = max((sv_ss_cid_T2470-sh_ss_cid_T2470));
[peak12,ind12] = min(abs((sv_ss_cid_T2470-sh_ss_cid_T2470)-peak12/2));
E50_ss_cid_T2470 = (sv_ss_cid_T2470(ind12)-sh_ss_cid_T2470(ind12))...
    /(ev_ss_cid_T2470(ind12));
pr_ss_cid_T2470 = (-eh_ss_cid_T2470(ind12))/(ev_ss_cid_T2470(ind12));

[peak13,~] = max((sv_ss_cid_T2472-sh_ss_cid_T2472));
[peak13,ind13] = min(abs((sv_ss_cid_T2472-sh_ss_cid_T2472)-peak13/2));
E50_ss_cid_T2472 = (sv_ss_cid_T2472(ind13)-sh_ss_cid_T2472(ind13))...
    /(ev_ss_cid_T2472(ind13));
pr_ss_cid_T2472 = (-eh_ss_cid_T2472(ind13))/(ev_ss_cid_T2472(ind13));

[peak14,~] = max((sv_ss_cid_T2476-sh_ss_cid_T2476));
[peak14,ind14] = min(abs((sv_ss_cid_T2476-sh_ss_cid_T2476)-peak14/2));
E50_ss_cid_T2476 = (sv_ss_cid_T2476(ind14)-sh_ss_cid_T2476(ind14))...
    /(ev_ss_cid_T2476(ind14));
pr_ss_cid_T2476 = (-eh_ss_cid_T2476(ind14))/(ev_ss_cid_T2476(ind14));

%% Relationship between dynamic and static E and PR from cores
% mean and std of dynamic measurements
% PRs
avg_pr_iso_shale_ciu_T2492 = mean(pr_iso_shale_ciu_T2492);
std_pr_iso_shale_ciu_T2492 = std(pr_iso_shale_ciu_T2492);
avg_pr_iso_shale_ciu_T2522 = mean(pr_iso_shale_ciu_T2522);
std_pr_iso_shale_ciu_T2522 = std(pr_iso_shale_ciu_T2522);

avg_pr_iso_ss_cid_T2461 = mean(pr_iso_ss_cid_T2461);
std_pr_iso_ss_cid_T2461 = std(pr_iso_ss_cid_T2461);
avg_pr_iso_ss_cid_T2463 = mean(pr_iso_ss_cid_T2463);
std_pr_iso_ss_cid_T2463 = std(pr_iso_ss_cid_T2463);
avg_pr_iso_ss_cid_T2464 = mean(pr_iso_ss_cid_T2464);
std_pr_iso_ss_cid_T2464 = std(pr_iso_ss_cid_T2464);
avg_pr_iso_ss_cid_T2466 = mean(pr_iso_ss_cid_T2466);
std_pr_iso_ss_cid_T2466 = std(pr_iso_ss_cid_T2466);
avg_pr_iso_ss_cid_T2467 = mean(pr_iso_ss_cid_T2467);
std_pr_iso_ss_cid_T2467 = std(pr_iso_ss_cid_T2467);
avg_pr_iso_ss_cid_T2469 = mean(pr_iso_ss_cid_T2469);
std_pr_iso_ss_cid_T2469 = std(pr_iso_ss_cid_T2469);
avg_pr_iso_ss_cid_T2470 = mean(pr_iso_ss_cid_T2470);
std_pr_iso_ss_cid_T2470 = std(pr_iso_ss_cid_T2470);
avg_pr_iso_ss_cid_T2472 = mean(pr_iso_ss_cid_T2472);
std_pr_iso_ss_cid_T2472 = std(pr_iso_ss_cid_T2472);
avg_pr_iso_ss_cid_T2476 = mean(pr_iso_ss_cid_T2476);
std_pr_iso_ss_cid_T2476 = std(pr_iso_ss_cid_T2476);

% Es
avg_E_iso_shale_ciu_T2492 = mean(E_iso_shale_ciu_T2492);
std_E_iso_shale_ciu_T2492 = std(E_iso_shale_ciu_T2492);
avg_E_iso_shale_ciu_T2522 = mean(E_iso_shale_ciu_T2522);
std_E_iso_shale_ciu_T2522 = std(E_iso_shale_ciu_T2522);

avg_E_iso_ss_cid_T2461 = mean(E_iso_ss_cid_T2461);
std_E_iso_ss_cid_T2461 = std(E_iso_ss_cid_T2461);
avg_E_iso_ss_cid_T2463 = mean(E_iso_ss_cid_T2463);
std_E_iso_ss_cid_T2463 = std(E_iso_ss_cid_T2463);
avg_E_iso_ss_cid_T2464 = mean(E_iso_ss_cid_T2464);
std_E_iso_ss_cid_T2464 = std(E_iso_ss_cid_T2464);
avg_E_iso_ss_cid_T2466 = mean(E_iso_ss_cid_T2466);
std_E_iso_ss_cid_T2466 = std(E_iso_ss_cid_T2466);
avg_E_iso_ss_cid_T2467 = mean(E_iso_ss_cid_T2467);
std_E_iso_ss_cid_T2467 = std(E_iso_ss_cid_T2467);
avg_E_iso_ss_cid_T2469 = mean(E_iso_ss_cid_T2469);
std_E_iso_ss_cid_T2469 = std(E_iso_ss_cid_T2469);
avg_E_iso_ss_cid_T2470 = mean(E_iso_ss_cid_T2470);
std_E_iso_ss_cid_T2470 = std(E_iso_ss_cid_T2470);
avg_E_iso_ss_cid_T2472 = mean(E_iso_ss_cid_T2472);
std_E_iso_ss_cid_T2472 = std(E_iso_ss_cid_T2472);
avg_E_iso_ss_cid_T2476 = mean(E_iso_ss_cid_T2476);
std_E_iso_ss_cid_T2476 = std(E_iso_ss_cid_T2476);

% fitting polynomial relationship for dynamic and static E modulus of sand
E_static_sand = [E50_ss_cid_T2461; E50_ss_cid_T2463; E50_ss_cid_T2464; E50_ss_cid_T2466;...
    E50_ss_cid_T2467; E50_ss_cid_T2469; E50_ss_cid_T2470; E50_ss_cid_T2472; E50_ss_cid_T2476];
E_dynamic_sand = [avg_E_iso_ss_cid_T2461; avg_E_iso_ss_cid_T2463; avg_E_iso_ss_cid_T2464; avg_E_iso_ss_cid_T2466;...
    avg_E_iso_ss_cid_T2467; avg_E_iso_ss_cid_T2469; avg_E_iso_ss_cid_T2470; avg_E_iso_ss_cid_T2472; avg_E_iso_ss_cid_T2476];

E_fit = polyfit(E_dynamic_sand,E_static_sand,1);
E_fit_val = polyval(E_fit,E_dynamic_sand);
E_equation = sprintf('E_{sta} = %.3f E_{dyn} - %.3f',E_fit(1),abs(E_fit(2)));

E_corr = corrcoef(E_static_sand,E_fit_val);
E_rsquared = E_corr(2)^2; E_r2 = num2str(E_rsquared);

%% SAVING ALL THE OUTPUTS
% save('rock_mechanics_data_outputs\all_data_v2.mat')