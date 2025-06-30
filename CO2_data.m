% Bulk Modulus of CO2 based on adiabatic eqns of PG-EoS
% Peng and Robinson, 1976 & Morse and Ingard, 1986 main papers/eqns
% Also see Carcione and Poletto, 2000 & Carcione and Picotti, 2006 & Picotti and Carcione, 2012
% Also see Batzle and Wang, 1992 & Wang and Nur, 1989
% Coded by Ufuk Durmus in 6/2024

% to see the PR-EoS model with Wang and Nur, 1989 data

clc; clear; close all;

%% INPUT PARAMETERS
% Digitized data of Wang and Nur, 1989
load('wang_nur_data\co2propdata.mat')

p = p./145.038; % pressure in MPa

%% PR-EoS
t_model = [27 57 97]; % temperature in C
p_model = 0:1:40;    % pressure in MPa

for i = 1:length(t_model)
    for j = 1:length(p_model)
[K_co2(i,j),rho_co2(i,j)] = CO2_CH4_rho_K_fun_final("CO2",t_model(i),p_model(j));
    end
end

Vp_co2 = (sqrt((K_co2)./(rho_co2))).*1000; % m/s

%% Save the results
% save('rock_physics_outputs\co2_outputs.mat')
%% RESULTS

figure
subplot(121)
plot(p_model,K_co2(:,:))
hold on
scatter(p,co2bulk(:,[2 5 9]),25,'k','filled')

subplot(122)
plot(p_model,rho_co2(:,:))
hold on
scatter(p,co2rho(:,[2 5 9]),25,'k','filled')