% Run script for AVA analysis in both isotropic and anisotropic media for PP and PS waves
% see also zoepp_edited, subZoeppritz, rueger_vti_avo_fun, aki_richards
% functions must be in the same directory
% coded by Ufuk Durmus in 12/2024
% updated with Graebner VTI functions for paper revision by Ufuk Durmus in 4/2025

close all; clear; clc;

%% Inputs/Models
p = 27.5; % pressure (MPa)
sw = 1; % water saturation
sgas = 1-sw; % gas saturation

p2 = 27.5; % pressure (MPa)
sw2 = 0.2; % water saturation
sgas2 = 1-sw2; % gas saturation

%% Estimate Fluid Properties for AVA Model
sal = 78400; % NaCl salinity (ppm)
t = 96; % temperature (C)

% brine properties
[K_brine,rho_brine] = batzle_wang_brine(sal,t,p);
[K_brine2,rho_brine2] = batzle_wang_brine(sal,t,p2);

% co2 properties
[K_co2,rho_co2] = CO2_CH4_rho_K_fun_final("CO2",t,p);
[K_co22,rho_co22] = CO2_CH4_rho_K_fun_final("CO2",t,p2);

% Fluid/Gas Mixture
rho_fluid = rho_brine.*sw + rho_co2.*sgas; % wood's mixing
rho_fluid2 = rho_brine2.*sw2 + rho_co22.*sgas2; % wood's mixing

% K_fluid_v = sgas.*K_co2 + sw.*K_brine; % voigt bound for patchy saturation

b_e = 5; % brie coefficient
K_fluid_brie = (K_brine-K_co2).*sw.^b_e+K_co2; % brie equation for patchy saturation
K_fluid_brie2 = (K_brine2-K_co22).*sw2.^b_e+K_co22; % brie equation for patchy saturation

%% Estimate the RP model for Sandstone AVA Model
% Johansen Fm. Mineralogy based on Sundal etal, 2016 (see the paper and XRD appendix excel file)
v_q = 0.67; v_feld = 0.18; v_plag = 0.03; v_clay = 0.12; % sum up to 1

phi = 0.26; % modeling porosity

% Use Voigt-Reuss-Hill bounds to estimate mineral mix
[K_voigt,K_reuss,K_vrh,mu_voigt,mu_reuss,mu_vrh,rho_matrix_rp,rho_dry_rp]= VRH_northernlights(v_q,v_feld,v_plag,v_clay,0,0,0,phi);

% Use Extended Maxwell Homogenization scheme with pore shape factor
s = 0.3; % pore shape factor

[K_eff,mu_eff]= Maxwell_iso_supersphere(K_vrh, mu_vrh, s, phi);

[K_sat,mu_sat,Vp_ss,Vs_ss,rhob] = gassmann_iso(K_eff,K_vrh,mu_eff,rho_dry_rp,phi,K_fluid_brie,rho_fluid);
[K_sat2,mu_sat2,Vp_ss2,Vs_ss2,rhob2] = gassmann_iso(K_eff,K_vrh,mu_eff,rho_dry_rp,phi,K_fluid_brie2,rho_fluid2);

%% Estimate the RP model for Shale AVA Model
% aspect ratios
ar_clay = 0.1;
ar_pore = 0.5;

% interclay medium
K_bw = 1.85; % assumption based on my previous work and Sayers papers 
G_bw = 0.4; % assumption based on my previous work and Sayers papers
v_bw = 0.08; % 8% for QEMSCAN to make it ~9% overall clay matrix

% shale properties
vclay_drake = 0.81; % volume of clay
phit_drake = 0.1; % shale porosity
rhob_sh = 2.5;

% mineralogy measurement type
min_type = 'QEMSCAN';

[C11_rp_ava,C33_rp_ava,C13_rp_ava,C55_rp_ava,C66_rp_ava,epsilon_rp_ava,gamma_rp_ava,delta_rp_ava,E1_rp_ava,E3_rp_ava,v13_rp_ava,v31_rp_ava,v12_rp_ava] = Maxwell_shale_well(K_bw,G_bw,ar_clay,ar_pore,vclay_drake,phit_drake,v_bw,min_type);

Vp_shale = (sqrt(C33_rp_ava./rhob_sh)).*1000; %m/s
Vs_shale = (sqrt(C55_rp_ava./rhob_sh)).*1000; %m/s

%% Define variables
vp1 = Vp_shale;
vs1 = Vs_shale;
rho1 = rhob_sh;
epsilon1 = epsilon_rp_ava;
gamma1 = gamma_rp_ava;
delta1 = delta_rp_ava;

vp2 = Vp_ss;
vs2 = Vs_ss;
rho2 = rhob;
epsilon2 = 0;
gamma2 = 0;
delta2 = 0;

vp22 = Vp_ss2;
vs22 = Vs_ss2;
rho22 = rhob2;

theta = 0:1:45;

%% ISOTROPIC EXACT ZOEPPRITZ EQUATIONS
[pp,ps,ss,sh,thetap,thetas]=zoepp_edited(vp1, vs1, rho1, vp2, vs2, rho2);
[pp2,ps2,ss2,sh2,thetap2,thetas2]=zoepp_edited(vp1, vs1, rho1, vp22, vs22, rho22);

%% VTI RUEGER RC APPROXIMATIONS
[R_vti_p,R_vti_ps,R_iso_ps,R_vti_sv,R_vti_sh]=rueger_vti_avo_fun(vp1,vs1,rho1,vp2,vs2,rho2,epsilon1,gamma1,delta1,epsilon2,gamma2,delta2,theta);
[R_vti_p2,R_vti_ps2,R_iso_ps2,R_vti_sv2,R_vti_sh2]=rueger_vti_avo_fun(vp1,vs1,rho1,vp22,vs22,rho22,epsilon1,gamma1,delta1,epsilon2,gamma2,delta2,theta);

%% VTI EXACT GRAEBNER EQUATIONS
[RP_vti, RPS_vti, thetap_vti] = Graebner_FUN_vel(vp1,vs1,rho1,epsilon1,delta1,vp2,vs2,rho2,epsilon2,delta2);
[RP_vti2, RPS_vti2, thetap_vti2] = Graebner_FUN_vel(vp1,vs1,rho1,epsilon1,delta1,vp22,vs22,rho22,epsilon2,delta2);

%% AKI & RICHARDS RC APPROXIMATION (just to see)
[RC] = aki_richards(vp1, vs1, rho1, vp2, vs2, rho2, theta);
[RC2] = aki_richards(vp1, vs1, rho1, vp22, vs22, rho22, theta);

%% SAVE THE NECESSARY OUTPUTS
% save('rock_physics_outputs\avo_rp_outputs.mat',"theta","thetap","thetap2",'pp','pp2'...
%     ,'R_vti_p','R_vti_p2','ps','ps2','R_vti_ps','R_vti_ps2','thetap_vti','thetap_vti2'...
%     ,'RP_vti','RP_vti2','RPS_vti','RPS_vti2')

%% PLOT THE RESULTS
figure
subplot(121)
plot(theta,zeros(length(theta)),'k','LineWidth',2)
hold on
plot(real(thetap),real(pp))
hold on
plot(real(thetap_vti),real(RP_vti),'-+')
hold on
plot(theta,R_vti_p,'--')
hold on
plot(theta,RC)
hold on
plot(real(thetap2),real(pp2))
hold on
plot(real(thetap_vti2),real(RP_vti2),'-+')
hold on
plot(theta,R_vti_p2,'--')
hold on
plot(theta,RC2)
xlim([0 40])
ylim([-0.3 0.101])

subplot(122)
plot(theta,zeros(length(theta)),'k','LineWidth',2)
hold on
plot(real(thetap),real(ps))
hold on
plot(real(thetap_vti),real(RPS_vti),'-+')
hold on
plot(theta,R_vti_ps,'--')
hold on
plot(theta,R_iso_ps)
hold on
plot(real(thetap_vti2),real(RPS_vti2),'-+')
hold on
plot(real(thetap2),real(ps2))
hold on
plot(theta,R_vti_ps2,'--')
hold on
plot(theta,R_iso_ps2)
xlim([0 40])
ylim([-0.3 0.101])