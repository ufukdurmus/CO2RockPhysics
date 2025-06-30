function [K_gas,rho_gas] = CO2_CH4_rho_K_fun_final(gas_type,t,p,v_co2,v_ch4)

% Thermodynamic properties of CO2, CH4 and air (PG-EoS)
% Based on Peng & Robinson, 1976 paper
% ANALYTICAL SOLUTION - ROOTS VERSION!!!
% Also see Carcione and Poletto, 2000 & Picotti and Carcione, 2012
% Coded by Ufuk Durmus in 6/2024

% INPUTS %
% t is temperature in C
% p is pressure in MPa
% gas_type index is the target gas (can be "CO2" "CH4" "air" or "CO2_CH4")
% v_co2 volume fraction of CO2 unitless (for mixing gas phase)
% v_ch4 volume fraction of CH4 unitless (for mixing gas phase)

% OUTPUT %
% rho is the density of gas at given temperature and pressure in g/cc
% K is the bulk modulus of gas at given temperature and pressure in GPa

%% DENSITY OF THE FLUID/GAS (CO2, CH4, air etc.)
ta = t + 273.15; % absolute temperature in K

% define chemical gas constants
r = 8.31; % gas constant J/mol K

if (gas_type == "CO2")
    
% for CO2
m_co2 = 44; % mass g/mole
w_co2 = 0.225; % acentric factor
tc_co2 = 31.1 + 273.15; % critical temperature in K
pc_co2 = 7.38; % critical pressure in MPa
% rhoc_co2 = 468.2; % critical density in kg/m3

% estimate Peng and Robinson Equation of State (PG-EoS)
% ta = t + 273.15; % absolute temperature in K
tr_co2 = ta./tc_co2; % reduced temperature unitless
pr_co2 = p./pc_co2; % reduced pressure unitless

beta_co2 = (1+(0.37464+1.54226.*w_co2-0.26992.*w_co2^2)*(1-sqrt(tr_co2))).^2;
a_co2 = (0.45724*r^2*tc_co2^2)/pc_co2*beta_co2;
b_co2 = (0.07780*r*tc_co2)/pc_co2;

A = a_co2*p/(r*ta).^2;
B = b_co2*p/(r*ta);
% Compressibility factor
Z = roots([1, -(1-B), (A-3.*B.^2-2.*B), -(A.*B-B.^2-B.^3)]);
% Choose the max Z for gas phase and min for liquid;
sol =(max((Z(imag(Z) == 0))));
% Calculate molar volume using Z = PV/RT
vm_co2 = (sol.*r.*ta)./p; % gas molar volume in cc/mol
rho_co2 = m_co2./vm_co2; % in g/cc

% heat capacity ratio at constant pressure (gamma)
% Piotti and Carcione modified this constant of Batzle and Wang, 1992 with
% Wang and Nur, 1989 data
gamma_co2 = 1.37 + (11.29./(pr_co2+6)) + (15.55./(pr_co2+1.3).^2)...
    - 38.89.*(exp(-1.25.*(pr_co2+1)));

% c = -(1/V)*(dV/dP) Gas Compressibility Eqn
derivative_co2 = -(r.*ta)./(vm_co2-b_co2).^2 +...
    ((2.*vm_co2+2.*b_co2).*a_co2)./(vm_co2.^2+2.*b_co2.*vm_co2-b_co2.^2).^2;

% K = 1/c = -V*(dV/dP)
Kt_co2 = -vm_co2.*derivative_co2; % isothermal bulk modulus
ct_co2 = 1./Kt_co2; % isothermal compressibility

% cs adiabatic compressibility & ct isothermal compressibility as cs=ct/gamma
cs_co2 = ct_co2./gamma_co2; % (MPa)^-1
K_co2 = 1./cs_co2; % MPa
K_co2 = K_co2./1000; % GPa

% results
rho_gas = rho_co2; % g/cc
K_gas = K_co2; % GPa

elseif (gas_type == "CH4")
    
% for CH4
m_ch4 = 16; % mass g/mole
w_ch4 = 0.0115; % acentric factor
tc_ch4 = -82.6 + 273.15; % critical temperature in K
pc_ch4 = 4.64; % critical pressure in MPa
% rhoc_ch4 = 162.7; % critical density in kg/m3
    
% estimate Peng and Robinson Equation of State (PG-EoS)
tr_ch4 = ta./tc_ch4; % reduced temperature unitless
pr_ch4 = p./pc_ch4; % reduced pressure

beta_ch4 = (1+(0.37464+1.54226.*w_ch4-0.26992.*w_ch4^2)*(1-sqrt(tr_ch4))).^2;
a_ch4 = (0.45724*r^2*tc_ch4^2)/pc_ch4*beta_ch4;
b_ch4 = (0.07780*r*tc_ch4)/pc_ch4;

A = a_ch4*p/(r*ta).^2;
B = b_ch4*p/(r*ta);
% Compressibility factor
Z = roots([1, -(1-B), (A-3*B^2-2*B), -(A*B-B^2-B^3)]);
% Choose the max Z for gas phase and min for liquid;
sol =(max((Z(imag(Z) == 0))));
% Calculate molar volume using Z = PV/RT
vm_ch4 = (sol*r*ta)/p; % gas molar volume in cc/mol
rho_ch4 = m_ch4./vm_ch4; % in g/cc

% heat capacity ratio at constant pressure (gamma) constant of Batzle and Wang, 1992
gamma_ch4 = 0.85+5.6./(pr_ch4+2)+27.1./(pr_ch4+3.5).^2-8.7.*exp(-0.65.*(pr_ch4+1));

% c = -(1/V)*(dV/dP) Gas Compressibility Eqn
derivative_ch4 = -(r.*ta)./(vm_ch4-b_ch4).^2 +...
    ((2.*vm_ch4+2.*b_ch4).*a_ch4)./(vm_ch4.^2+2.*b_ch4.*vm_ch4-b_ch4.^2).^2;

% K = 1/c = -V*(dV/dP)
Kt_ch4 = -vm_ch4.*derivative_ch4; % isothermal bulk modulus
ct_ch4 = 1./Kt_ch4; % isothermal compressibility

% cs adiabatic compressibility & ct isothermal compressibility as cs=ct/gamma
cs_ch4 = ct_ch4./gamma_ch4; % (MPa)^-1
K_ch4 = 1./cs_ch4; % MPa
K_ch4 = K_ch4./1000; % GPa

% results
rho_gas = rho_ch4; % g/cc
K_gas = K_ch4; % GPa

elseif (gas_type == "air")

% for air
m_air = 29; % mass g/mole
w_air = 0.078; % acentric factor
tc_air = -140.8 + 273.15; % critical temperature in K
pc_air = 3.7; % critical pressure in MPa
% rhoc_air = 340; % critical density in kg/m3

% estimate Peng and Robinson Equation of State (PG-EoS)
tr_air = ta./tc_air; % reduced temperature unitless
%pr_air = ta./pc_air; % reduced pressure

beta_air = (1+(0.37464+1.54226.*w_air-0.26992.*w_air^2)*(1-sqrt(tr_air))).^2;
a_air = (0.45724*r^2*tc_air^2)/pc_air*beta_air;
b_air = (0.07780*r*tc_air)/pc_air;

A = a_air*p/(r*ta).^2;
B = b_air*p/(r*ta);
% Compressibility factor
Z = roots([1, -(1-B), (A-3*B^2-2*B), -(A*B-B^2-B^3)]);
% Choose the max Z for gas phase and min for liquid;
sol =(max((Z(imag(Z) == 0))));
% Calculate molar volume using Z = PV/RT
vm_air = (sol*r*ta)/p; % gas molar volume in cc/mol
rho_air = m_air./vm_air; % in g/cc

% heat capacity ratio at constant pressure (gamma) constant
gamma_air = 1.4; % many references 1.4 for dry air, can be 4/3 from Carcione

% c = -(1/V)*(dV/dP) Gas Compressibility Eqn
derivative_air = -(r.*ta)./(vm_air-b_air).^2 +...
    ((2.*vm_air+2.*b_air).*a_air)./(vm_air.^2+2.*b_air.*vm_air-b_air.^2).^2;

% K = 1/c = -V*(dV/dP)
Kt_air = -vm_air.*derivative_air; % isothermal bulk modulus
ct_air = 1./Kt_air; % isothermal compressibility

% cs adiabatic compressibility & ct isothermal compressibility as cs=ct/gamma
cs_air = ct_air./gamma_air; % (MPa)^-1
K_air = 1./cs_air; % MPa
K_air = K_air./1000; % GPa

% results
rho_gas = rho_air; % g/cc
K_gas = K_air; % GPa

elseif (nargin>3 && gas_type == "CO2_CH4")
    
m_mix = [44 16]; % mass g/mole
w_mix = [0.225 0.0115]; % acentric factor
tc_mix= [31.1 + 273.15  -82.6 + 273.15]; % critical temperature in K
pc_mix = [7.38  4.64]; % critical pressure in MPa
v_avg = [v_co2 v_ch4] ; % volume fractions
avg_m = sum(m_mix.*v_avg); % avg molar mass
tr_mix = ta./tc_mix;

pr_co2 = p./pc_mix(1); % reduced pressure
pr_ch4 = p./pc_mix(2); % reduced pressure

beta_co2_ch4 = (1+(0.37464+1.54226.*w_mix-0.26992.*w_mix.^2).*(1-sqrt(tr_mix))).^2;
a = 0.45724.*((r.*tc_mix).^2./pc_mix).*beta_co2_ch4;
b = 0.0778*r.*tc_mix./pc_mix;
[a_co2_ch4,b_co2_ch4] = VanDerWaals_Mixing(a,b,v_avg);

A = a_co2_ch4*p/(r*ta).^2;
B = b_co2_ch4*p/(r*ta);
% Compressibility factor
Z = roots([1, -(1-B), (A-3*B^2-2*B), -(A*B-B^2-B^3)]);
% Choose the max Z for gas phase and min for liquid;
sol =(max((Z(imag(Z) == 0))));
% Calculate molar volume using Z = PV/RT
vm_co2_ch4 = (sol*r*ta)/p; % gas molar volume in cc/mol
rho_co2_ch4 = avg_m./vm_co2_ch4; % in g/cc

% heat capacity ratio at constant pressure (gamma) constant
% picotti and carcione, 2012 & Danesh, 2001 PR-EoS mixing rules
gamma_co2 = 1.37 + (11.29./(pr_co2+6)) + (15.55./(pr_co2+1.3).^2)...
    - 38.89.*(exp(-1.25.*(pr_co2+1)));

gamma_ch4 = 0.85+5.6./(pr_ch4+2)+27.1./(pr_ch4+3.5).^2-8.7.*exp(-0.65.*(pr_ch4+1));
% mixing
gamma_co2_ch4 = gamma_co2.*(v_co2.^2) + gamma_ch4.*(v_ch4.^2) + ...
    v_co2.*v_ch4.*(gamma_co2+gamma_ch4);

% c = -(1/V)*(dV/dP) Gas Compressibility Eqn
derivative_co2_ch4 = -(r.*ta)./(vm_co2_ch4-b_co2_ch4).^2 +...
    ((2.*vm_co2_ch4+2.*b_co2_ch4).*a_co2_ch4)./(vm_co2_ch4.^2+2.*b_co2_ch4.*vm_co2_ch4-b_co2_ch4.^2).^2;

% K = 1/c = -V*(dV/dP)
Kt_co2_ch4 = -vm_co2_ch4.*derivative_co2_ch4; % isothermal bulk modulus
ct_co2_ch4 = 1./Kt_co2_ch4; % isothermal compressibility

% cs adiabatic compressibility & ct isothermal compressibility as cs=ct/gamma
cs_co2_ch4 = ct_co2_ch4./gamma_co2_ch4; % (MPa)^-1
K_co2_ch4 = 1./cs_co2_ch4; % MPa
K_co2_ch4 = K_co2_ch4./1000; % GPa

% results
rho_gas = rho_co2_ch4; % g/cc
K_gas = K_co2_ch4; % GPa

elseif(nargin <=3 && gas_type == "CO2_CH4")
    error('Please insert volume fractions of CO2 and CH4 in the mixture!')

else
    error('gas type must be CO2 CH4 air or CO2_CH4!')
end
end