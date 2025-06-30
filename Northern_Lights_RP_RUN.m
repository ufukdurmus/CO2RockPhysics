%% ROCK PHYSICS ANALYSIS OF NORTHERN LIGHTS FOR CO2 STORAGE
% see also rock_mechanics_data, well_logs_data
% VRH_northernlights, Maxwell_shale, Maxwell_iso_supersphere 
% should be in the same directory
% coded by Ufuk Durmus in 12/2024

close all; clear; clc;

%% LOAD THE NECESSARY DATA
% inputs from processed rock_mechanics_data
load('rock_mechanics_data_outputs\all_data_v2.mat')

% inputs from processed well logs data
load('well_data_outputs\well_rp_outputs.mat')

%% SS ISOTROPIC ROCK PHYSICS MODEL
% Johansen Fm. Mineralogy based on Sundal etal, 2016 (see the paper and XRD appendix excel file)
v_q = 0.67; v_feld = 0.18; v_plag = 0.03; v_clay = 0.12; % sum up to 1

% TEST the model with T2464 sample Johansen-2 Fm.
v_por_ss = phi_initial(12); % porosity of T2464 sample

% Use Voigt-Reuss-Hill bounds to estimate mineral mix
[K_voigt,K_reuss,K_vrh,mu_voigt,mu_reuss,mu_vrh,rho_matrix_rp,rho_dry_rp]= VRH_northernlights(v_q,v_feld,v_plag,v_clay,0,0,0,v_por_ss);

% Use Extended Maxwell Homogenization scheme with pore shape factor
s = 0.25:0.05:0.5; % pore shape factor
s_por = 0:0.05:0.4; % modeling porosity

for j = 1:length(s)
    for i = 1:length(s_por)
[K_eff(i,j),mu_eff(i,j)]= Maxwell_iso_supersphere(K_vrh, mu_vrh, s(j), s_por(i));
    end
end

% plots; all the data with rp model
figure
subplot(1,2,1)
scatter(phi_initial(10),c33_ss_cid_T2461,'k','filled')
hold on
scatter(phi_initial(13),c33_ss_cid_T2466,'b','filled')
hold on
scatter(phi_initial(16),c33_ss_cid_T2470,'r','filled')
hold on
scatter(phi_initial(11),c33_ss_cid_T2463,'k','filled')
hold on
scatter(phi_initial(14),c33_ss_cid_T2467,'b','filled')
hold on
scatter(phi_initial(17),c33_ss_cid_T2472,'r','filled')
hold on
scatter(phi_initial(12),c33_ss_cid_T2464,'k','filled')
hold on
scatter(phi_initial(15),c33_ss_cid_T2469,'b','filled')
hold on
scatter(phi_initial(18),c33_ss_cid_T2476,'r','filled')
hold on
scatter(phi_initial(19),c33_ss_hydro_T2479,'k','filled')
hold on
scatter(phi_initial(20),c33_ss_hydro_T2480,'b','filled')
hold on
scatter(phi_initial(21),c33_ss_hydro_T2482,'r','filled')
hold on
scatter(phi_initial(7),c33_ss_caust_T2505,'k','filled')
hold on
% scatter(phi_initial(8),c33_ss_caust_T2507,'b','filled') % no axial data
hold on
scatter(phi_initial(9),c33_ss_caust_T2511,'r','filled')
hold on
plot(s_por,K_eff+(4/3).*mu_eff,'--k')
xlabel('Porosity','FontSize',12,'FontWeight','bold')
ylabel('P-wave Modulus (GPa)','FontSize',12,'FontWeight','bold')
ylim([0 inf])
% xlim([10 35])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';

subplot(1,2,2)
pl1 = scatter(phi_initial(10),c55_ss_cid_T2461,'k','filled');
hold on
pl2 = scatter(phi_initial(13),c55_ss_cid_T2466,'b','filled');
hold on
pl3 = scatter(phi_initial(16),c55_ss_cid_T2470,'r','filled');
hold on
scatter(phi_initial(11),c55_ss_cid_T2463,'k','filled')
hold on
scatter(phi_initial(14),c55_ss_cid_T2467,'b','filled')
hold on
scatter(phi_initial(17),c55_ss_cid_T2472,'r','filled')
hold on
scatter(phi_initial(12),c55_ss_cid_T2464,'k','filled')
hold on
scatter(phi_initial(15),c55_ss_cid_T2469,'b','filled')
hold on
scatter(phi_initial(18),c55_ss_cid_T2476,'r','filled')
hold on
scatter(phi_initial(19),c55_ss_hydro_T2479,'k','filled')
hold on
scatter(phi_initial(20),c55_ss_hydro_T2480,'b','filled')
hold on
scatter(phi_initial(21),c55_ss_hydro_T2482,'r','filled')
hold on
scatter(phi_initial(7),c55_ss_caust_T2505,'k','filled')
hold on
scatter(phi_initial(8),c55_ss_caust_T2507,'b','filled')
hold on
scatter(phi_initial(9),c55_ss_caust_T2511,'r','filled')
hold on
plot(s_por,mu_eff,'--k')
xlabel('Porosity','FontSize',12,'FontWeight','bold')
ylabel('S-wave Modulus (GPa)','FontSize',12,'FontWeight','bold')
ylim([0 inf])
% xlim([10 35])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
legend([pl1(1),pl2(1),pl3(1)],'Johansen-2','Johansen-4','Cook-2','Fontweight','bold')

%% SHALE ANISOTROPIC ROCK PHYSICS MODEL

% ROCK PHYSICS MODELING OF THE WELL LOGS
% for anisotropic shale
% aspect ratios
ar_clay = 0.1;
ar_pore = 0.5;

% interclay medium
K_bw1 = 1.85; % assumption based on my previous work and Sayers papers 
G_bw1 = 0.4; % assumption based on my previous work and Sayers papers
v_bw1 = 0.08; % 8% for QEMSCAN and 6% for XRD mineralogy to make it ~9% overall clay matrix

K_bw2 = 1.85;
G_bw2 = 0.4;
v_bw2 = 0.06;

K_bw3 = 1.85;
G_bw3 = 0.4;
v_bw3 = 0.08; % 0.0909 when normalized with clay from Intra-Drake for info not for modeling

K_bw4 = 1.55;
G_bw4 = 0.4;
v_bw4 = 0.035;

K_bw5 = 1.55;
G_bw5 = 0.4;
v_bw5 = 0.041;

K_bw6 = 1.55;
G_bw6 = 0.4;
v_bw6 = 0.047; % 0.0629 when normalized with clay from Lower-Drake1 for info not for modeling

%%
v_interparticle = v_bw1; % makes it about 9 percent
v_si = 0.08;  % volume fraction of Inclusion 1 (1n our case "Smectite")
v_i = 0.62;  % volume fraction of Inclusion 2 (1n our case "Illite")
v_k = 0.1;  % volume fraction of Inclusion 3 (1n our case "Kaolinite")
v_c = 0;  % volume fraction of Inclusion 4 (1n our case "Chlorite")

sum_clay = v_interparticle+v_si+v_i+v_k+v_c;

% normalized clay mineral composition
v_si = v_si/sum_clay;  
v_i = v_i/sum_clay; 
v_k = v_k/sum_clay;
v_c = v_c/sum_clay;
v_interparticle = v_interparticle./sum_clay;

%%

% mineralogy measurement type
min_type = 'QEMSCAN';
min_type_lowdrake1 = 'QEMSCAN_lowdrake1';

% Drake shale Fm is divided into two parts for modeling purposes

for t = 1:length(vclay_drake)
[C11_rp1(t),C33_rp1(t),C13_rp1(t),C55_rp1(t),C66_rp1(t),epsilon_rp1(t),gamma_rp1(t),delta_rp1(t),E1_rp1(t),E3_rp1(t),v13_rp1(t),v31_rp1(t),v12_rp1(t)] = Maxwell_shale_well(K_bw1,G_bw1,ar_clay,ar_pore,vclay_drake(t),phit_drake(t),v_bw1,min_type);
end

for t2 = 1:length(vclay_drake2)
[C11_rp2(t2),C33_rp2(t2),C13_rp2(t2),C55_rp2(t2),C66_rp2(t2),epsilon_rp2(t2),gamma_rp2(t2),delta_rp2(t2),E1_rp2(t2),E3_rp2(t2),v13_rp2(t2),v31_rp2(t2),v12_rp2(t2)] = Maxwell_shale_well(K_bw2,G_bw2,ar_clay,ar_pore,vclay_drake2(t2),phit_drake2(t2),v_bw2,min_type);
end

for t3 = 1:length(vclay_drake3)
[C11_rp3(t3),C33_rp3(t3),C13_rp3(t3),C55_rp3(t3),C66_rp3(t3),epsilon_rp3(t3),gamma_rp3(t3),delta_rp3(t3),E1_rp3(t3),E3_rp3(t3),v13_rp3(t3),v31_rp3(t3),v12_rp3(t3)] = Maxwell_shale_well(K_bw3,G_bw3,ar_clay,ar_pore,vclay_drake3(t3),phit_drake3(t3),v_bw3,min_type);
end

for t4 = 1:length(vclay_lowdrake1)
[C11_rp4(t4),C33_rp4(t4),C13_rp4(t4),C55_rp4(t4),C66_rp4(t4),epsilon_rp4(t4),gamma_rp4(t4),delta_rp4(t4),E1_rp4(t4),E3_rp4(t4),v13_rp4(t4),v31_rp4(t4),v12_rp4(t4)] = Maxwell_shale_well(K_bw4,G_bw4,ar_clay,ar_pore,vclay_lowdrake1(t4),phit_lowdrake1(t4),v_bw4,min_type_lowdrake1);
end

for t5 = 1:length(vclay_lowdrake12)
[C11_rp5(t5),C33_rp5(t5),C13_rp5(t5),C55_rp5(t5),C66_rp5(t5),epsilon_rp5(t5),gamma_rp5(t5),delta_rp5(t5),E1_rp5(t5),E3_rp5(t5),v13_rp5(t5),v31_rp5(t5),v12_rp5(t5)] = Maxwell_shale_well(K_bw5,G_bw5,ar_clay,ar_pore,vclay_lowdrake12(t5),phit_lowdrake12(t5),v_bw5,min_type_lowdrake1);
end

for t6 = 1:length(vclay_lowdrake13)
[C11_rp6(t6),C33_rp6(t6),C13_rp6(t6),C55_rp6(t6),C66_rp6(t6),epsilon_rp6(t6),gamma_rp6(t6),delta_rp6(t6),E1_rp6(t6),E3_rp6(t6),v13_rp6(t6),v31_rp6(t6),v12_rp6(t6)] = Maxwell_shale_well(K_bw6,G_bw6,ar_clay,ar_pore,vclay_lowdrake13(t6),phit_lowdrake13(t6),v_bw6,min_type_lowdrake1);
end

% merging the shale parts after modeling;
% rp model results merged
C11_rp_shale = [C11_rp1 C11_rp2 C11_rp3 C11_rp4 C11_rp5 C11_rp6];
C33_rp_shale = [C33_rp1 C33_rp2 C33_rp3 C33_rp4 C33_rp5 C33_rp6];
C13_rp_shale = [C13_rp1 C13_rp2 C13_rp3 C13_rp4 C13_rp5 C13_rp6];
C55_rp_shale = [C55_rp1 C55_rp2 C55_rp3 C55_rp4 C55_rp5 C55_rp6];
C66_rp_shale = [C66_rp1 C66_rp2 C66_rp3 C66_rp4 C66_rp5 C66_rp6];
epsilon_rp_shale = [epsilon_rp1 epsilon_rp2 epsilon_rp3 epsilon_rp4 epsilon_rp5 epsilon_rp6];
gamma_rp_shale = [gamma_rp1 gamma_rp2 gamma_rp3 gamma_rp4 gamma_rp5 gamma_rp6];
delta_rp_shale = [delta_rp1 delta_rp2 delta_rp3 delta_rp4 delta_rp5 delta_rp6];
E1_rp_shale = [E1_rp1 E1_rp2 E1_rp3 E1_rp4 E1_rp5 E1_rp6];
E3_rp_shale = [E3_rp1 E3_rp2 E3_rp3 E3_rp4 E3_rp5 E3_rp6];
v13_rp_shale = [v13_rp1 v13_rp2 v13_rp3 v13_rp4 v13_rp5 v13_rp6];
v31_rp_shale = [v31_rp1 v31_rp2 v31_rp3 v31_rp4 v31_rp5 v31_rp6];
v12_rp_shale = [v12_rp1 v12_rp2 v12_rp3 v12_rp4 v12_rp5 v12_rp6];
% well log data merged
C33_well_shale = [C33_drake; C33_drake2; C33_drake3; C33_lowdrake1; C33_lowdrake12; C33_lowdrake13];
C55_well_shale = [C55_drake; C55_drake2; C55_drake3; C55_lowdrake1; C55_lowdrake12; C55_lowdrake13];
depth_well_shale = [depth_drake; depth_drake2; depth_drake3; depth_lowdrake1; depth_lowdrake12; depth_lowdrake13];
phit_well_shale = [phit_drake; phit_drake2; phit_drake3; phit_lowdrake1; phit_lowdrake12; phit_lowdrake13];
vclay_well_shale = [vclay_drake; vclay_drake2; vclay_drake3; vclay_lowdrake1; vclay_lowdrake12; vclay_lowdrake13];

figure
subplot(1,5,1)
plot(C11_rp_shale,depth_well_shale,'k','LineWidth',2)
hold on
errorbar(mean(c11_shale_ciu_T2529),2592.75,std(c11_shale_ciu_T2529),'horizontal',"o","MarkerSize",2,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',2)
set(gca,'YDir','reverse')

subplot(1,5,2)
plot(C33_rp_shale,depth_well_shale,'k','LineWidth',2)
hold on
plot(C33_well_shale,depth_well_shale,'r','LineWidth',2)
hold on
errorbar(mean(c33_shale_ciu_T2522),2592.65,std(c33_shale_ciu_T2522),'horizontal',"o","MarkerSize",2,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',2)
set(gca,'YDir','reverse')
xlim([10 40])

subplot(1,5,3)
plot(C55_rp_shale,depth_well_shale,'k','LineWidth',2)
hold on
plot(C55_well_shale,depth_well_shale,'r','LineWidth',2)
hold on
errorbar(mean(c55_shale_ciu_T2492),2592.81,std(c55_shale_ciu_T2492),'horizontal',"o","MarkerSize",2,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',2)
set(gca,'YDir','reverse')
xlim([3 12])

subplot(1,5,4)
plot(C66_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

subplot(1,5,5)
plot(C13_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

figure
subplot(1,4,1)
plot(epsilon_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

subplot(1,4,2)
plot(gamma_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

subplot(1,4,3)
plot(delta_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

subplot(1,4,4)
plot(vclay_well_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

figure
subplot(1,5,1)
plot(E1_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

subplot(1,5,2)
plot(E3_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

subplot(1,5,3)
plot(v13_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

subplot(1,5,4)
plot(v31_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

subplot(1,5,5)
plot(v12_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')

%% for isotropic sandstone
% Johansen Fm. and Cook Fm. Mineralogy based on EoS well logs
vq_sand = 1-vclay_sand; % sum up to 1

% Use Voigt-Reuss-Hill bounds to estimate mineral mix
[K_v_sand,K_r_sand,K_vrh_sand,mu_v_sand,mu_r_sand,mu_vrh_sand]= VRH_northernlights(vq_sand,0,0,vclay_sand,0,0,0,0);

% Use Extended Maxwell Homogenization scheme with pore shape factor
s_sand = 0.3; % pore shape factor

for z = 1:length(vclay_sand)
[K_sand(z),mu_sand(z)]= Maxwell_iso_supersphere(K_vrh_sand(z), mu_vrh_sand(z), s_sand, phit_sand(z));
end
C33_rp_sand = (K_sand+(4/3).*mu_sand);

% Gassmann Fluid Substitution
K_fluid = 2.85;
rho_fluid = 1.03;

for z = 1:length(vclay_sand)
[K_sat(z)] = gassmann_iso(K_sand(z),K_vrh_sand(z),mu_sand(z),[],phit_sand(z),K_fluid,rho_fluid);
end

C33_sat_sand = (K_sat+(4/3).*mu_sand);
C33_sat_sand(isnan(C33_sat_sand)) = 60; % some values are NaN in the array

% values are spiky for both well log and model results
C33_sat_sand(C33_sat_sand>42) = 42;
C33_sand(C33_sand>42) = 42;

mu_sand(mu_sand>15) = 15;
C55_sand(C55_sand>15) = 15;

% FOR RMS ERROR CALCULATIONS OF SANDSTONE INTERVALS
C33_sat_sand=C33_sat_sand'; C33_sat_sand([349:395 697:743])=C33_sand([349:395 697:743]);
mu_sand=mu_sand'; mu_sand([349:395 697:743])=C55_sand([349:395 697:743]);

figure
subplot(1,2,1)
plot(C33_sat_sand,depth_sand,'k','LineWidth',2)
hold on
plot(C33_sand,depth_sand,'r','LineWidth',2)
set(gca,'YDir','reverse')
xlim([10 50])

subplot(1,2,2)
plot(mu_sand,depth_sand,'k','LineWidth',2)
hold on
plot(C55_sand,depth_sand,'r','LineWidth',2)
set(gca,'YDir','reverse')
xlim([0 20])

%% SAVE THE NECESSARY OUTPUTS
% save('rock_physics_outputs\rp_outputs.mat',"C11_rp_shale","C33_rp_shale","C13_rp_shale","C55_rp_shale","C66_rp_shale"...
%     ,"C33_sat_sand","mu_sand","epsilon_rp_shale","gamma_rp_shale","delta_rp_shale","E1_rp_shale","E3_rp_shale"...
%     ,"v13_rp_shale","v31_rp_shale","v12_rp_shale","C33_well_shale","C55_well_shale","depth_well_shale"...
%     ,"s_por","K_eff","mu_eff")