% Run script for 4D reservoir response analysis
% functions must be in the same directory
% coded by Ufuk Durmus in 12/2024
% updated for paper revision by Ufuk Durmus in 4/2025

close all; clear; clc;

%% Inputs/Models
p = 27.5; % pressure (MPa)
p2 = 37.5; % pressure (MPa)
sal = 78400; % NaCl salinity (ppm)
t = 96; % temperature (C)

phi = 0.0001:0.001:0.3501; % modeling porosity
sw = 1:-0.001:0; % water saturation
sgas = 1-sw; % gas saturation

for i = 1:length(phi)
    for j = 1:length(sw)

%% Estimate Fluid Properties for AVA Model

% brine properties
[K_brine,rho_brine] = batzle_wang_brine(sal,t,p);

% co2 properties
[K_co2,rho_co2] = CO2_CH4_rho_K_fun_final("CO2",t,p);

% Fluid/Gas Mixture
rho_fluid(:,j) = rho_brine.*sw(j) + rho_co2.*sgas(j); % wood's mixing
% K_fluid_v(:,j) = sgas(j).*K_co2 + sw(j).*K_brine; % voigt bound for patchy saturation

b_e = 5; % brie coefficient
K_fluid_brie(:,j) = (K_brine-K_co2).*sw(j).^b_e+K_co2; % brie equation for patchy saturation

%% Estimate the RP model for Sandstone AVA Model
% Johansen Fm. Mineralogy based on Sundal etal, 2016 (see the paper and XRD appendix excel file)
v_q = 0.67; v_feld = 0.18; v_plag = 0.03; v_clay = 0.12; % sum up to 1
% v_q = 1; v_feld = 0; v_plag = 0; v_clay = 0; % pure sandstone deneme

% Use Voigt-Reuss-Hill bounds to estimate mineral mix
[K_voigt(i),K_reuss(i),K_vrh(i),mu_voigt(i),mu_reuss(i),mu_vrh(i),rho_matrix_rp(i),rho_dry_rp(i)]= VRH_northernlights(v_q,v_feld,v_plag,v_clay,0,0,0,phi(i));

% Use Extended Maxwell Homogenization scheme with pore shape factor
s = 0.3; % pore shape factor

[K_eff(i),mu_eff(i)]= Maxwell_iso_supersphere(K_vrh(i), mu_vrh(i), s, phi(i));

[K_sat(i,j),mu_sat(i,j),Vp_ss(i,j),Vs_ss(i,j),rhob(i,j)] = gassmann_iso(K_eff(i),K_vrh(i),mu_eff(i),rho_dry_rp(i),phi(i),K_fluid_brie(j),rho_fluid(j));

    end
end

%% 2ND PART-CASE
for i = 1:length(phi)
    for j = 1:length(sw)

%% Estimate Fluid Properties for AVA Model

% brine properties
[K_brine2,rho_brine2] = batzle_wang_brine(sal,t,p2);

% co2 properties
[K_co22,rho_co22] = CO2_CH4_rho_K_fun_final("CO2",t,p2);

% Fluid/Gas Mixture
rho_fluid2(:,j) = rho_brine2.*sw(j) + rho_co22.*sgas(j); % wood's mixing
% K_fluid_v2(:,j) = sgas(j).*K_co22 + sw(j).*K_brine2; % voigt bound for patchy saturation

b_e = 5; % brie coefficient
K_fluid_brie2(:,j) = (K_brine2-K_co22).*sw(j).^b_e+K_co22; % brie equation for patchy saturation

%% Estimate the RP model for Sandstone AVA Model
% Johansen Fm. Mineralogy based on Sundal etal, 2016 (see the paper and XRD appendix excel file)
v_q = 0.67; v_feld = 0.18; v_plag = 0.03; v_clay = 0.12; % sum up to 1
% v_q = 1; v_feld = 0; v_plag = 0; v_clay = 0; % pure sandstone deneme

% Use Voigt-Reuss-Hill bounds to estimate mineral mix
[K_voigt2(i),K_reuss2(i),K_vrh2(i),mu_voigt2(i),mu_reuss2(i),mu_vrh2(i),rho_matrix_rp2(i),rho_dry_rp2(i)]= VRH_northernlights(v_q,v_feld,v_plag,v_clay,0,0,0,phi(i));

% Use Extended Maxwell Homogenization scheme with pore shape factor
s = 0.3; % pore shape factor

[K_eff2(i),mu_eff2(i)]= Maxwell_iso_supersphere(K_vrh2(i), mu_vrh2(i), s, phi(i));

[K_sat2(i,j),mu_sat2(i,j),Vp_ss2(i,j),Vs_ss2(i,j),rhob2(i,j)] = gassmann_iso(K_eff2(i),K_vrh2(i),mu_eff2(i),rho_dry_rp2(i),phi(i),K_fluid_brie2(j),rho_fluid2(j));

    end
end
%%
Zp = Vp_ss.*rhob;
Zs = Vs_ss.*rhob;
vpvs = Zp./Zs;

Zp2 = Vp_ss2.*rhob2;
Zs2 = Vs_ss2.*rhob2;

deltazp = Zp-Zp2;
deltazs = Zs-Zs2;

%% Save
% save('rock_physics_outputs\imp_phi_co2_rp_outputs.mat',"Zp","Zs","sgas","phi")

%% Plots

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
