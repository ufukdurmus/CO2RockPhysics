function [K_voigt,K_reuss,K_vrh,mu_voigt,mu_reuss,mu_vrh,rho_matrix_rp,rho_dry_rp]= VRH_northernlights(v_q,v_feld,v_plag,v_clay,v_cal,v_dol,v_pyr,v_por)

% Voigt Reuss Hill Average Bounds

% Coded by Ufuk Durmus
% Last Updated 3/6/2019

% Rock Physics Handbook by Gary Mavko et al, page 110-111

%% Define Mineralogy and Volume Fractions

% Carmichael, R. S., 1989, Practical handbook of physical properties of
% rocks and minerals: CRC Press.

K_q = 37;
m_q = 44;

K_cal = 70.2;
m_cal = 29;

K_dol = 94.9;
m_dol = 45;

K_feld = 37.5;
m_feld = 15;

K_plag = 75.6;
m_plag = 25.6;

% K_sid = 123.7;
% m_sid = 51;

K_pyrite = 147.4;
m_pyrite = 132.5;

% Vanorio T., Prasad M. and Nur A. 2003. Elastic properties of dry
% clay mineral aggregates, suspensions and sandstones. Geophysical
% Journal International 155, 319–326.

% K_clay = 12;
% m_clay = 6;

% K_clay = 21;
% m_clay = 7;

K_clay = 25;
m_clay = 9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcuting Voigt-Reuss-Hill Elastic Moduli (isotropic) average

K_voigt = ((v_q.*K_q)+(v_feld.*K_feld)+(v_plag.*K_plag)+(v_cal.*K_cal)...
    +(v_pyr.*K_pyrite)+(v_dol.*K_dol)+(v_clay.*K_clay));

K_reuss = 1./((v_q./K_q)+(v_feld./K_feld)+(v_plag./K_plag)+(v_cal./K_cal)...
    +(v_pyr./K_pyrite)+(v_dol./K_dol)+(v_clay./K_clay));

K_vrh = (K_voigt+K_reuss)./2;

mu_voigt = ((v_q.*m_q)+(v_feld.*m_feld)+(v_plag.*m_plag)+(v_cal.*m_cal)...
    +(v_pyr.*m_pyrite)+(v_dol.*m_dol)+(v_clay.*m_clay));

mu_reuss = 1./((v_q./m_q)+(v_feld./m_feld)+(v_plag./m_plag)+(v_cal./m_cal)...
    +(v_pyr./m_pyrite)+(v_dol./m_dol)+(v_clay./m_clay));

mu_vrh = (mu_voigt+mu_reuss)./2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional dry density estimation for QC with core measurements

% Carmichael, R. S., 1989, Practical handbook of physical properties of
% rocks and minerals: CRC Press.

rho_q = 2.65;

rho_cal = 2.71;

rho_dol = 2.87;

rho_feld = 2.62;

rho_plag = 2.63;

% rho_sid = 3.96;

rho_pyrite = 4.93;

% Vanorio T., Prasad M. and Nur A. 2003. Elastic properties of dry
% clay mineral aggregates, suspensions and sandstones. Geophysical
% Journal International 155, 319–326.

rho_clay = 2.55;

% rho_clay = 2.60;

% DENSITY ESTIMATION

rho_matrix_rp = (v_q.*rho_q)+(v_feld.*rho_feld)+(v_plag.*rho_plag)+(v_cal.*rho_cal)...
    +(v_pyr.*rho_pyrite)+(v_dol.*rho_dol)+(v_clay.*rho_clay);

rho_dry_rp = rho_matrix_rp.*(1-v_por);

end