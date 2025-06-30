% function to model well log data with Extended Maxwell Scheme
% for Intra-Drake shale in Northern Lights

function [C11,C33,C13,C55,C66,epsilon,gamma,delta,E_1,E_3,v13,v31,v12] = Maxwell_shale_well(K_interparticle,G_interparticle,gamma_clay,gamma_por,v_clay,v_por,v_bw,clay_type)

%% Inclusion Properties - Constant Inputs

gamma1 = 0.99;    % aspect ratio of Inclusion 1 (for spheroids = 100 and for d1sks = 0.01)
gamma2 = gamma_por;    % aspect ratio of Inclusion 2 (for spheroids = 100 and for d1sks = 0.01)
gamma3 = 0.99;    % aspect ratio of Inclusion 3 (for spheroids = 100 and for d1sks = 0.01)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clay Properties - Constant Inputs

% aspect ratios
gamma_si = gamma_clay;    % aspect ratio of Inclusion 1 (for spheroids = 100 and for d1sks = 0.01)
gamma_i = gamma_clay;    % aspect ratio of Inclusion 2 (for spheroids = 100 and for d1sks = 0.01)
gamma_k = gamma_clay;    % aspect ratio of Inclusion 3 (for spheroids = 100 and for d1sks = 0.01)
gamma_c = gamma_clay;    % aspect ratio of Inclusion 4 (for spheroids = 100 and for d1sks = 0.01)

% bulk and shear moduli of clays
K_si = 12.3;
G_si = 15.6;

K_i = 60.2;
G_i = 25.4;

K_k = 55.5;
G_k = 31.8;

K_c = 54.3;
G_c = 30.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MINERALOGY for Intra-Drake shale from Grande etal, 2024,Thomson etal, 2024 ARMA PAPERS (3 of them), Rahman etal, 2022
% there are 2 methods reported for mineralogy

% XRD Analysis
v_q = 1-v_clay;                    % volume fraction of quartz
v_q(v_q<0) = 0;
v_cal = 0;                         % calcite not needed for this study

sum_min = (v_q+v_cal+v_clay+v_por); 

% Normalize each component by dividing by the sum of minerals

v_q = v_q./(sum_min);                    
v_cal = v_cal./(sum_min);                  
v_clay = v_clay./(sum_min);                
v_por = v_por./(sum_min);

% Inclusions
if (clay_type == "XRD") % background Cijs from Maxwell method directly estimated from XRD mineralogy
% from XRD analysis
v_interparticle = v_bw; % assumption based on my previous work and Sayers papers, makes it about 9 percent
v_si = 0;  % volume fraction of Inclusion 1 (1n our case "Smectite")
v_i = 0.44;  % volume fraction of Inclusion 2 (1n our case "Illite")
v_k = 0;  % volume fraction of Inclusion 3 (1n our case "Kaolinite")
v_c = 0.17;  % volume fraction of Inclusion 4 (1n our case "Chlorite")

sum_clay = v_interparticle+v_si+v_i+v_k+v_c;

% normalized clay mineral composition
v_si = v_si/sum_clay;  
v_i = v_i/sum_clay; 
v_k = v_k/sum_clay;
v_c = v_c/sum_clay;
v_interparticle = v_interparticle./sum_clay;

% Background Estimation (Clay Matrix)
    
[C0_11,C0_33,C0_13,C0_55,C0_66]= Maxwell_clay (K_interparticle,G_interparticle,K_i,G_i,K_k,G_k,K_c,G_c,v_i,v_k,v_c,gamma_i,gamma_k,gamma_c);

elseif (clay_type == "QEMSCAN") % background Cijs from Maxwell method directly estimated from QEMSCAN mineralogy
% from QEMSCAN analysis
v_interparticle = v_bw; % assumption based on my previous work and Sayers papers, makes it about 9 percent
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

% Background Estimation (Clay Matrix)
    
[C0_11,C0_33,C0_13,C0_55,C0_66]= Maxwell_clay (K_interparticle,G_interparticle,K_i,G_i,K_k,G_k,K_si,G_si,v_i,v_k,v_si,gamma_i,gamma_k,gamma_c);

elseif (clay_type == "QEMSCAN_lowdrake1") % background Cijs from Maxwell method directly estimated from QEMSCAN mineralogy
% from QEMSCAN analysis
v_interparticle = v_bw; % assumption based on my previous work and Sayers papers, makes it about 9 percent
v_si = 0.05;  % volume fraction of Inclusion 1 (1n our case "Smectite")
v_i = 0.45;  % volume fraction of Inclusion 2 (1n our case "Illite")
v_k = 0.2;  % volume fraction of Inclusion 3 (1n our case "Kaolinite")
v_c = 0;  % volume fraction of Inclusion 4 (1n our case "Chlorite")

sum_clay = v_interparticle+v_si+v_i+v_k+v_c;

% normalized clay mineral composition
v_si = v_si/sum_clay;  
v_i = v_i/sum_clay; 
v_k = v_k/sum_clay;
v_c = v_c/sum_clay;
v_interparticle = v_interparticle./sum_clay;

% Background Estimation (Clay Matrix)
    
[C0_11,C0_33,C0_13,C0_55,C0_66]= Maxwell_clay (K_interparticle,G_interparticle,K_i,G_i,K_k,G_k,K_si,G_si,v_i,v_k,v_si,gamma_i,gamma_k,gamma_c);


elseif (clay_type == "literature") % background Cijs from literature

% illite mica aggragate by Vernik and Kachanov, 2010 & Vernik, 2016(book)
% C0_11 = 53.4;
% C0_33 = 33.4;
% C0_13 = 21;
% C0_55 = 8.5;
% C0_66 = 12.7;

% Ortega, 2007
C0_11 = 44.9;
C0_33 = 24.2;
C0_13 = 18.1;
C0_55 = 3.7;
C0_66 = 11.6;

% Keller, 2017
% C0_11 = 24.77;
% C0_33 = 15.86;
% C0_13 = 7.7;
% C0_55 = 5.26;
% C0_66 = 7.47;

% trial from core measurements of rock mechanics data
% C0_11 = 43;
% % C0_11 = 46.5;
% C0_33 = 21.5;
% % C0_33 = 26.5;
% C0_13 = 18;
% C0_55 = 5;
% C0_66 = 8;
% % C0_55 = 7;
% % C0_66 = 12;

else

    error('clay background estimation should be XRD or QEMSCAN!')

end

%% Computation of properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maxwell Homogenization Scheme
K1 = 37;
G1 = 44;

K2 = 2.85; % from Batzle and Wang with ~PVT properties
G2 = 0;

K3 = 70.2;
G3 = 29;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inclusions & Aspect ratios
v1 = v_q;  % volume fraction of Inclusion 1 (1n our case "quartz")
v2 = v_por;  % volume fraction of Inclusion 2 (1n our case "fluid-filled pores")
v3 = v_cal;  % volume fraction of Inclusion 3 (1n our case "calcite") % not needed for this project

%% RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[C11,C33,C13,C55,C66]= Maxwell_ani_shale(C0_11,C0_33,C0_13,C0_55,C0_66,K1,G1,K2,G2,K3,G3,v1,v2,v3,gamma1,gamma2,gamma3);
C12 = C11-2.*C66;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Thomsen Parameters

epsilon = (C11-C33)/(2.*C33);
gamma = (C66-C55)/(2.*C55);
delta = ((C13+C55).^2-(C33-C55).^2)./((2.*C33).*(C33-C55));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating Young's Modulus and Poisson's Ratio in Anisotropic Media

E_3 = (C33.*(C11+C12)-2.*(C13^2))./(C11+C12);

E_1 = C11 + (C13.^2.*(-2.*C11+C12)+C12.*(-C33.*C12+C13.^2))...
./(C33.*C11-C13.^2);

v31 = C13./(C11+C12);

v12 = (C33.*C12-C13.^2)./(C11.*C33-C13.^2);

v13 = (C13.*(C11-C12))./(C11.*C33-C13.^2);

end

