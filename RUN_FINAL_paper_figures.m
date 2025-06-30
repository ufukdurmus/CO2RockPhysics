% ROCK PHYSICS MODELING OF 31/5-7 EOS WELL FOR CO2 STORAGE
% This script is used for creating the final figures for the paper after extensive analysis and modeling
% also see rock_mechanics_data_v2, well_logs_data, Northern_Lights_RP_RUN, CO2_data
% Coded by Ufuk Durmus on 12/2024

close all; clear; clc;

%% LOAD THE NECESSARY DATA
% inputs from processed rock_mechanics_data
load('rock_mechanics_data_outputs\all_data_v2.mat')

% inputs from processed well logs data
load('well_data_outputs\well_rp_outputs.mat')

% inputs from rock physics modeling
load('rock_physics_outputs\rp_outputs.mat')

% inputs from rock physics AVA modeling
load('rock_physics_outputs\avo_rp_outputs.mat')

% inputs from rock physics 4D modeling
load('rock_physics_outputs\imp_phi_co2_rp_outputs.mat')

% inputs from CO2 property modeling
load('rock_physics_outputs\co2_outputs.mat')

%% FIGURE 1: Showing well logs
ylimit = [2580 2880];
well_logs_disp_final(depth,depth2,depth_core_well,ylimit,cali,gr,vp,vs,rhob,rhob_core_well,drho,phit,phit_core_well,v_clay,intra_drake,lower_drake1,cook4,cook2,burton,johansen4,johansen2,amundsen)

%% FIGURE : Velocity vs Pressure from Ultrasonic Core Measurements
figure
subplot(121)
scatter(sdiff_shale_ciu_T2492,vp_axial_shale_ciu_T2492,20,"black",'filled','o')
hold on
scatter(sdiff_shale_ciu_T2492,vs_axial_shale_ciu_T2492,20,"blue",'filled','o')
hold on
scatter(sdiff_shale_ciu_T2514,vp_axial_shale_ciu_T2514,20,"black",'filled','o')
hold on
scatter(sdiff_shale_ciu_T2522,vp_axial_shale_ciu_T2522,20,"black",'filled','o')
hold on
scatter(sdiff_shale_ciu_T2522,vs_axial_shale_ciu_T2522,20,"blue",'filled','o')
hold on
scatter(sdiff_shale_ciu_T2529,vp_axial_shale_ciu_T2529,20,"red",'filled','o')
hold on
scatter(sdiff_shale_ciu_T2533,vp_axial_shale_ciu_T2533,20,"red",'filled','o')
xlabel('\sigma_{diff} (\sigma_v^,-\sigma_h^,) (MPa)','FontSize',12,'FontWeight','bold')
ylabel('Velocity (m/s)','FontSize',12,'FontWeight','bold')
ylim([1000 4500])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
legend('V_{P0}','V_{S0}','','','','V_{P90}','Location','northwest','Fontweight','bold','Location','east')
annotation("textbox",[.03 .8 .1 .2],'String','a)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(122)
scatter(sdiff_ss_cid_T2464,vp_axial_ss_cid_T2464,20,'k','filled','o')
hold on
scatter(sdiff_ss_cid_T2464,vs_axial_ss_cid_T2464,20,'k','filled','^')
hold on
scatter(sdiff_ss_cid_T2469,vp_axial_ss_cid_T2469,20,'b','filled','o')
hold on
scatter(sdiff_ss_cid_T2469,vs_axial_ss_cid_T2469,20,'b','filled','^')
hold on
scatter(sdiff_ss_cid_T2472,vp_axial_ss_cid_T2472,20,'r','filled','o')
hold on
scatter(sdiff_ss_cid_T2472,vs_axial_ss_cid_T2472,20,'r','filled','^')
xlabel('\sigma_{diff} (\sigma_v^,-\sigma_h^,) (MPa)','FontSize',12,'FontWeight','bold')
ylabel('Velocity (m/s)','FontSize',12,'FontWeight','bold')
ylim([1000 4500])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
legend('Johansen-2','','Johansen-4','','Cook-2','','Location','northwest','Fontweight','bold','Location','east')
annotation("textbox",[.47 .8 .1 .2],'String','b)','FontWeight','bold','FontSize',14,'EdgeColor','None')

%% FIGURE : CO2 properties modeling with data
figure
subplot(121)
plot(p_model,K_co2(:,:),'LineWidth',2)
hold on
scatter(p,co2bulk(:,2),30,'k','o','filled')
hold on
scatter(p,co2bulk(:,5),30,'k','s','filled')
hold on
scatter(p,co2bulk(:,9),30,'k','^','filled')
xlabel('Pressure (MPa)','FontSize',12,'FontWeight','bold')
ylabel('K_{CO_2} (GPa)','FontSize',12,'FontWeight','bold')
ylim([0 inf])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
legend('27 ^\circC','57 ^\circC','97 ^\circC','Location','northwest','Fontweight','bold')
annotation("textbox",[.05 .8 .1 .2],'String','a)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(122)
plot(p_model,rho_co2(:,:),'LineWidth',2)
hold on
scatter(p,co2rho(:,2),30,'k','o','filled')
hold on
scatter(p,co2rho(:,5),30,'k','s','filled')
hold on
scatter(p,co2rho(:,9),30,'k','^','filled')
xlabel('Pressure (MPa)','FontSize',12,'FontWeight','bold')
ylabel('\rho_{CO_2} (g/cc)','FontSize',12,'FontWeight','bold')
ylim([0 inf])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
annotation("textbox",[.5 .8 .1 .2],'String','b)','FontWeight','bold','FontSize',14,'EdgeColor','None')

%% FIGURE : RP model calibration for sandstone using core
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
ylim([0 80])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
annotation("textbox",[.04 .80 .1 .2],'String','a)','FontWeight','bold','FontSize',16,'EdgeColor','None')

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
ylim([0 30])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
legend([pl1(1),pl2(1),pl3(1)],'Johansen-2','Johansen-4','Cook-2','Fontweight','bold')
annotation("textbox",[.48 .80 .1 .2],'String','b)','FontWeight','bold','FontSize',16,'EdgeColor','None')

%% FIGURE : RP model calibration for sandstone using well
figure
subplot(1,2,1)
plot(C33_sat_sand,depth_sand,'k','LineWidth',2)
hold on
plot(C33_sand,depth_sand,'r','LineWidth',2)
set(gca,'YDir','reverse')
xlabel('M (GPa)','FontSize',12,'FontWeight','bold')
ylabel('Depth (m)','FontSize',12,'FontWeight','bold')
grid on
xlim([10 50])
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.05 .81 .1 .2],'String','a)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,2,2)
plot(mu_sand,depth_sand,'k','LineWidth',2)
hold on
plot(C55_sand,depth_sand,'r','LineWidth',2)
set(gca,'YDir','reverse')
xlabel('\mu (GPa)','FontSize',12,'FontWeight','bold')
set(gca,'Yticklabel',[])
grid on
xlim([2 18])
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.50 .81 .1 .2],'String','b)','FontWeight','bold','FontSize',14,'EdgeColor','None')

% visualize the tops
r_tops = 2:2:18; % showing the tops in the range as in vclay section
hold on
plot(r_tops,cook2.*ones(size(r_tops)),'black','LineWidth',2)
text(18,cook2,'cook2','Color','black','FontWeight','bold')
hold on
plot(r_tops,burton.*ones(size(r_tops)),'black','LineWidth',2)
text(18,burton,'burton','Color','black','FontWeight','bold')
hold on
plot(r_tops,johansen4.*ones(size(r_tops)),'black','LineWidth',2)
text(18,johansen4,'johansen4','Color','black','FontWeight','bold')
hold on
plot(r_tops,johansen2.*ones(size(r_tops)),'black','LineWidth',2)
text(18,johansen2,'johansen2','Color','black','FontWeight','bold')

%% FIGURE : RP model calibration for shale
figure
subplot(1,8,1)
plot(C11_rp_shale,depth_well_shale,'k','LineWidth',2)
hold on
errorbar(mean(c11_shale_ciu_T2529),2592.75,std(c11_shale_ciu_T2529),'horizontal',"o","MarkerSize",2,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',2)
set(gca,'YDir','reverse')
xlabel('C_{11} (GPa)','FontSize',12,'FontWeight','bold')
ylabel('Depth (m)','FontSize',12,'FontWeight','bold')
grid on
xlim([35 65])
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.10 .81 .1 .2],'String','a)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,8,2)
plot(C33_rp_shale,depth_well_shale,'k','LineWidth',2)
hold on
plot(C33_well_shale,depth_well_shale,'r','LineWidth',2)
hold on
errorbar(mean(c33_shale_ciu_T2522),2592.65,std(c33_shale_ciu_T2522),'horizontal',"o","MarkerSize",2,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlim([10 40])
xlabel('C_{33} (GPa)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.21 .81 .1 .2],'String','b)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,8,3)
plot(C55_rp_shale,depth_well_shale,'k','LineWidth',2)
hold on
plot(C55_well_shale,depth_well_shale,'r','LineWidth',2)
hold on
errorbar(mean(c55_shale_ciu_T2492),2592.81,std(c55_shale_ciu_T2492),'horizontal',"o","MarkerSize",2,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlim([2 12])
xlabel('C_{55} (GPa)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.31 .81 .1 .2],'String','c)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,8,4)
plot(C66_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlim([12 24])
xlabel('C_{66} (GPa)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.41 .81 .1 .2],'String','d)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,8,5)
plot(C13_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlim([8 14])
xlabel('C_{13} (GPa)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.51 .81 .1 .2],'String','e)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,8,6)
plot(epsilon_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlabel('\epsilon','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.61 .81 .1 .2],'String','f)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,8,7)
plot(gamma_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlabel('\gamma','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.70 .81 .1 .2],'String','g)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,8,8)
plot(delta_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlabel('\delta','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.80 .81 .1 .2],'String','h)','FontWeight','bold','FontSize',14,'EdgeColor','None')

% visualize the tops
r_tops = -0.08:0.001:-0.06; % showing the tops in the range as in vclay section
hold on
plot(r_tops,(intra_drake+5).*ones(size(r_tops)),'black','LineWidth',2)
text(-0.06,intra_drake+5,'intra-drake','Color','black','FontWeight','bold')
hold on
plot(r_tops,lower_drake1.*ones(size(r_tops)),'black','LineWidth',2)
text(-0.06,lower_drake1,'lower-drake1','Color','black','FontWeight','bold')

%% FIGURE : Stress vs Strain
% SHALES
figure
subplot(121)
plot(ev_shale_ciu_T2492,sv_shale_ciu_T2492-sh_shale_ciu_T2492,'k','LineWidth',2)
hold on
plot(ev_shale_ciu_T2514,sv_shale_ciu_T2514-sh_shale_ciu_T2514,'k','LineWidth',2)
hold on
plot(ev_shale_ciu_T2522,sv_shale_ciu_T2522-sh_shale_ciu_T2522,'k','LineWidth',2)
hold on
plot(ev_shale_ciu_T2529,sv_shale_ciu_T2529-sh_shale_ciu_T2529,'r','LineWidth',2)
hold on
plot(ev_shale_ciu_T2533,sv_shale_ciu_T2533-sh_shale_ciu_T2533,'r','LineWidth',2)
hold on
plot(eh_shale_ciu_T2492,sv_shale_ciu_T2492-sh_shale_ciu_T2492,'k--','LineWidth',2)
hold on
plot(eh_shale_ciu_T2514,sv_shale_ciu_T2514-sh_shale_ciu_T2514,'k--','LineWidth',2)
hold on
plot(eh_shale_ciu_T2522,sv_shale_ciu_T2522-sh_shale_ciu_T2522,'k--','LineWidth',2)
hold on
plot(eh_shale_ciu_T2529,sv_shale_ciu_T2529-sh_shale_ciu_T2529,'r--','LineWidth',2)
hold on
plot(eh_shale_ciu_T2533,sv_shale_ciu_T2533-sh_shale_ciu_T2533,'r--','LineWidth',2)
hold on
xlabel('Strain \epsilon_h and  \epsilon_v (mS)','FontSize',12,'FontWeight','bold')
ylabel('Differential Stress, \sigma_{diff} (MPa)','FontSize',12,'FontWeight','bold')
ylim([0 45])
xlim([-25 25])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
legend('vertical','','', 'horizontal', '','FontWeight','bold','Location','northwest')
annotation("textbox",[.09 .78 .1 .2],'String','a)','FontWeight','bold','FontSize',16,'EdgeColor','None')

% SANDSTONES

subplot(122)
plot(ev_ss_cid_T2461,sv_ss_cid_T2461-sh_ss_cid_T2461,'g','LineWidth',2)
hold on
plot(ev_ss_cid_T2463,sv_ss_cid_T2463-sh_ss_cid_T2463,'Color','[0.4660 0.6740 0.1880]','LineWidth',2)
hold on
plot(ev_ss_cid_T2464,sv_ss_cid_T2464-sh_ss_cid_T2464,'Color','[0 0.5 0]','LineWidth',2)
hold on
plot(ev_ss_cid_T2466,sv_ss_cid_T2466-sh_ss_cid_T2466,'r','LineWidth',2)
hold on
plot(ev_ss_cid_T2467,sv_ss_cid_T2467-sh_ss_cid_T2467,'Color','[0.8500 0.3250 0.0980]','LineWidth',2)
hold on
plot(ev_ss_cid_T2469,sv_ss_cid_T2469-sh_ss_cid_T2469,'Color','[0.6350 0.0780 0.1840]','LineWidth',2)
hold on
plot(ev_ss_cid_T2470,sv_ss_cid_T2470-sh_ss_cid_T2470,'Color','[0.3010 0.7450 0.9330]','LineWidth',2)
hold on
plot(ev_ss_cid_T2472,sv_ss_cid_T2472-sh_ss_cid_T2472,'Color','[0 0.4470 0.7410]','LineWidth',2)
hold on
plot(ev_ss_cid_T2476,sv_ss_cid_T2476-sh_ss_cid_T2476,'b','LineWidth',2)
hold on
plot(eh_ss_cid_T2461,sv_ss_cid_T2461-sh_ss_cid_T2461,'g','LineWidth',2,'LineStyle','--')
hold on
plot(eh_ss_cid_T2463,sv_ss_cid_T2463-sh_ss_cid_T2463,'Color','[0.4660 0.6740 0.1880]','LineWidth',2,'LineStyle','--')
hold on
plot(eh_ss_cid_T2464,sv_ss_cid_T2464-sh_ss_cid_T2464,'Color','[[0 0.5 0]]','LineWidth',2,'LineStyle','--')
hold on
plot(eh_ss_cid_T2466,sv_ss_cid_T2466-sh_ss_cid_T2466,'r--','LineWidth',2)
hold on
plot(eh_ss_cid_T2467,sv_ss_cid_T2467-sh_ss_cid_T2467,'Color','[0.8500 0.3250 0.0980]','LineWidth',2,'LineStyle','--')
hold on
plot(eh_ss_cid_T2469,sv_ss_cid_T2469-sh_ss_cid_T2469,'Color','[0.6350 0.0780 0.1840]','LineWidth',2,'LineStyle','--')
hold on
plot(eh_ss_cid_T2470,sv_ss_cid_T2470-sh_ss_cid_T2470,'Color','[0.3010 0.7450 0.9330]','LineWidth',2,'LineStyle','--')
hold on
plot(eh_ss_cid_T2472,sv_ss_cid_T2472-sh_ss_cid_T2472,'Color','[0 0.4470 0.7410]','LineWidth',2,'LineStyle','--')
hold on
plot(eh_ss_cid_T2476,sv_ss_cid_T2476-sh_ss_cid_T2476,'b','LineWidth',2,'LineStyle','--')
xlabel('Strain \epsilon_h and  \epsilon_v (mS)','FontSize',12,'FontWeight','bold')
ylabel('Differential Stress, \sigma_{diff} (MPa)','FontSize',12,'FontWeight','bold')
ylim([0 100])
xlim([-50 50])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
legend('','','Johansen-2', 'Johansen-4', '','','','','Cook-2','FontWeight','bold','Location','northwest')
annotation("textbox",[.52 .78 .1 .2],'String','b)','FontWeight','bold','FontSize',16,'EdgeColor','None')

%% FIGURE : E static vs dynamic
figure
errorbar(avg_E_iso_ss_cid_T2461,E50_ss_cid_T2461,std_E_iso_ss_cid_T2461,'horizontal',"o","MarkerSize",10,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',0.75)
hold on
errorbar(avg_E_iso_ss_cid_T2463,E50_ss_cid_T2463,std_E_iso_ss_cid_T2463,'horizontal',"o","MarkerSize",10,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',0.75)
hold on
errorbar(avg_E_iso_ss_cid_T2464,E50_ss_cid_T2464,std_E_iso_ss_cid_T2464,'horizontal',"o","MarkerSize",10,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',0.75)
hold on
errorbar(avg_E_iso_ss_cid_T2466,E50_ss_cid_T2466,std_E_iso_ss_cid_T2466,'horizontal',"^","MarkerSize",10,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',0.75)
hold on
errorbar(avg_E_iso_ss_cid_T2467,E50_ss_cid_T2467,std_E_iso_ss_cid_T2467,'horizontal',"^","MarkerSize",10,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',0.75)
hold on
errorbar(avg_E_iso_ss_cid_T2469,E50_ss_cid_T2469,std_E_iso_ss_cid_T2469,'horizontal',"^","MarkerSize",10,...
    "Color",'b',"MarkerFaceColor",'b','LineWidth',0.75)
hold on
errorbar(avg_E_iso_ss_cid_T2470,E50_ss_cid_T2470,std_E_iso_ss_cid_T2470,'horizontal',"s","MarkerSize",10,...
    "Color",'r',"MarkerFaceColor",'r','LineWidth',0.75)
hold on
errorbar(avg_E_iso_ss_cid_T2472,E50_ss_cid_T2472,std_E_iso_ss_cid_T2472,'horizontal',"s","MarkerSize",10,...
    "Color",'r',"MarkerFaceColor",'r','LineWidth',0.75)
hold on
errorbar(avg_E_iso_ss_cid_T2476,E50_ss_cid_T2476,std_E_iso_ss_cid_T2476,'horizontal',"s","MarkerSize",10,...
    "Color",'r',"MarkerFaceColor",'r','LineWidth',0.75)
hold on
plot(E_dynamic_sand,E_fit_val,'k','LineWidth',2)
hold on
text(11,19,E_equation,'FontWeight','bold')
hold on
text(11,17,['R^2 =' sprintf('% .3f ',E_rsquared)],'FontWeight','bold')
xlabel('E_{dynamic} (GPa)','FontSize',12,'FontWeight','bold')
ylabel('E_{static} (GPa)','FontSize',12,'FontWeight','bold')
ylim([0 20])
xlim([10 35])
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
legend('Johansen-2','','','Johansen-4','','','Cook-2','Location','southeast','FontWeight','bold')

%% FIGURE : Geomechanics plot; E, PR, K0
% anisotropic shale
figure
subplot(1,4,1)
plot(E1_rp_shale,depth_well_shale,'k','LineWidth',2)
set(gca,'YDir','reverse')
xlabel('E_1 (GPa)','FontSize',12,'FontWeight','bold')
ylabel('Depth (m)','FontSize',12,'FontWeight','bold')
xlim([25 55])
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.07 .81 .1 .2],'String','a)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,4,2)
plot(E3_rp_shale,depth_well_shale,'k','LineWidth',2)
hold on
errorbar(avg_E_iso_shale_ciu_T2522,2592.65,std_E_iso_shale_ciu_T2522,'horizontal',"o","MarkerSize",3,...
    "Color",'r',"MarkerFaceColor",'r','LineWidth',0.75)
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlabel('E_3 (GPa)','FontSize',12,'FontWeight','bold')
xlim([15 35])
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.29 .81 .1 .2],'String','b)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,4,3)
plot(v13_rp_shale,depth_well_shale,'k','LineWidth',2)
hold on
scatter(pr_shale_ciu_T2529,2592.75,25,'b','filled')
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlabel('\nu_{13}','FontSize',12,'FontWeight','bold')
xlim([0.25 0.4])
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.5 .81 .1 .2],'String','c)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(1,4,4)
plot(v31_rp_shale,depth_well_shale,'k','LineWidth',2)
hold on
errorbar(avg_pr_iso_shale_ciu_T2522,2592.65,std_pr_iso_shale_ciu_T2522,'horizontal',"o","MarkerSize",3,...
    "Color",'r',"MarkerFaceColor",'r','LineWidth',0.75)
hold on
scatter(pr_shale_ciu_T2514,2592.65,25,'b','filled')
hold on
scatter(pr_shale_ciu_T2522,2592.65,25,'b','filled')
set(gca,'YDir','reverse')
set(gca,'Yticklabel',[])
xlabel('\nu_{31}','FontSize',12,'FontWeight','bold')
xlim([0.14 0.22])
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.XAxisLocation = 'top';
ax.XGrid = 'on';
ax.YGrid = 'off';
annotation("textbox",[.7 .81 .1 .2],'String','d)','FontWeight','bold','FontSize',14,'EdgeColor','None')

%% FIGURE : AVA analysis iso vs aniso
figure
subplot(121)
plot(theta,zeros(length(theta)),'k','LineWidth',1)
hold on
plot(real(thetap),real(pp),'b','LineWidth',2)
hold on
% plot(theta,R_vti_p,'b--','LineWidth',2) % previous Rueger results
plot(real(thetap_vti),real(RP_vti),'b--','LineWidth',2)
hold on
plot(real(thetap2),real(pp2),'r','LineWidth',2)
hold on
% plot(theta,R_vti_p2,'r--','LineWidth',2) % previous Rueger results
plot(real(thetap_vti2),real(RP_vti2),'r--','LineWidth',2)
xlim([0 40])
ylim([-0.301 0.051])
ylabel('P-wave reflection coefficient','FontSize',12,'FontWeight','bold')
xlabel('angle of incidence (^\circ)','FontSize',12,'FontWeight','bold')
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
annotation("textbox",[.03 .8 .1 .2],'String','a)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(122)
plot(theta,zeros(length(theta)),'k','LineWidth',1);
hold on
pt1 = plot(real(thetap),real(ps),'b','LineWidth',2);
hold on
% pt2 = plot(theta,R_vti_ps,'b--','LineWidth',2); % previous Rueger results
pt2 = plot(real(thetap_vti),real(RPS_vti),'b--','LineWidth',2);
hold on
pt3 = plot(real(thetap2),real(ps2),'r','LineWidth',2);
hold on
% pt4 = plot(theta,R_vti_ps2,'r--','LineWidth',2); % previous Rueger results
pt4 = plot(real(thetap_vti2),real(RPS_vti2),'r--','LineWidth',2);
xlim([0 40])
ylim([-0.301 0.051])
ylabel('PS-wave reflection coefficient','FontSize',12,'FontWeight','bold')
xlabel('angle of incidence (^\circ)','FontSize',12,'FontWeight','bold')
grid on
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
legend([pt1(1),pt2(1),pt3(1),pt4(1)],'base-iso','base-vti','monitor-iso', 'monitor-vti','Fontweight','bold','Location','southeast')
annotation("textbox",[.48 .8 .1 .2],'String','b)','FontWeight','bold','FontSize',14,'EdgeColor','None')

%% FIGURE : 4D plot PHIT and Sco2

figure
subplot(121)
imagesc(sgas,phi,Zp./1000)
xlabel('S_{CO2} (v/v)','FontSize',12,'FontWeight','bold')
ylabel('Porosity (v/v)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.YDir = 'normal';
c=colorbar;
c.Label.String = ('Zp (km/s*g/cc)');
c.FontWeight = ('bold');
c.FontSize = (12);
cmap = (jet);
colormap(cmap)
clim([5 13])
annotation("textbox",[.05 .8 .1 .2],'String','a)','FontWeight','bold','FontSize',14,'EdgeColor','None')

subplot(122)
imagesc(sgas,phi,Zs./1000)
xlabel('S_{CO2} (v/v)','FontSize',12,'FontWeight','bold')
%ylabel('Porosity (v/v)','FontSize',12,'FontWeight','bold')
ax = gca;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontWeight = 'bold';
ax.YDir = 'normal';
c=colorbar;
c.Label.String = ('Zs (km/s*g/cc)');
c.FontWeight = ('bold');
c.FontSize = (12);
cmap = (jet);
colormap(cmap)
clim([3 8])
annotation("textbox",[.5 .8 .1 .2],'String','b)','FontWeight','bold','FontSize',14,'EdgeColor','None')
