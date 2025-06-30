% Gassmann Isotropic Fluid Substitution Code
% coded by Ufuk Durmus on 02/27/2019, revisited on 7/2024
% based on Mavko, 2009 and Smith, 2003

%% FUNCTION
function [K_sat,mu_sat,Vp_sat,Vs_sat,rhob] = gassmann_iso(K_eff,K_matrix,mu_eff,rho_eff,phi,K_fluid,rho_fluid)
% note that I name K_frame in the paper as  K_eff (the rock with porosity)

%% Saturated Rock Density - Bulk Density
rhob = rho_eff + rho_fluid.*phi; %rho_eff = rho_frame, rock matrix with porosity = rho0*(1-phi)

%% Fluid Substitution for Saturated Bulk Modulus
K_sat = K_eff + (((1-(K_eff./K_matrix)).^2)./((phi./K_fluid))+((1-phi)./K_matrix)-(K_eff./(K_matrix.^2)));

%% Saturated Vp and Vs values
mu_sat = mu_eff;

Vp_sat = (sqrt((K_sat+(4/3).*mu_sat)./rhob)).*1000; %m/s
Vs_sat = (sqrt(mu_sat./rhob)).*1000; %m/s

end
