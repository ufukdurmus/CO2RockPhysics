function [K_brine,rho_brine]=batzle_wang_brine(sal,t,p)
% Batzle-Wang relations only for brine properties

% INPUTS %
% sal: NaCl Sal (ppm)
% p: pore pressure (Mpa)
% t: Rock temperature (C)

% Kreuss: Reuss bound of mixed fluid's bulk modulus, homogeneous saturation
% Kvoigt: Voigt bount of mixed fluid's bulk modulus, patchy saturation

% OUTPUTS %
% rho_brine: Density of brine in g/cm3
% K_brine: Bulk Moduli of brine in Gpa

% Rearranged by Ufuk Durmus in 7/2024

%% Input parameters
% salinity is in ppm divided by 1e6
sal = sal./1e6;

%% brine density===========================================================
rhow = 1+1e-6.*(-80.*t-3.3.*t.^2+0.00175.*t.^3+489.*p-2.*t.*p+ ...
0.016.*t.^2.*p-1.3e-5.*t.^3.*p-0.333.*p.^2-0.002.*t.*p.^2);

rho_brine = rhow+sal.*(0.668+0.44.*sal+1e-6.*(300.*p-2400.*p.*sal+ ...
t.*(80+3.*t-3300.*sal-13.*p+47.*p.*sal))); %g/cc

%% brine velocity==========================================================
% water coeefficient matrix
matrixw = [1402.85 4.871 -0.04783 1.487e-4 -2.197e-7
1.524 -0.0111 2.747e-4 -6.503e-7 7.987e-10
3.437e-3 1.739e-4 -2.135e-6 -1.455e-8 5.230e-11
-1.197e-5 -1.628e-6 1.237e-8 1.327e-10 -4.614e-13]';

% water velocity
velw=0;
for i=1:5
for j=1:4
  velw=velw+matrixw(i,j).*t.^(i-1).*p.^(j-1);
end
end

% gas-free brine Bulk Modulus
vpb0 = velw+sal.*(1170-9.6.*t+0.055.*t.^2-8.5e-5.*t.^3+2.6.*p- ...
0.0029.*t.*p-0.0476.*p.^2)+sal.^1.5.*(780-10.*p+0.16.*p.^2)-1820.*sal.^2;

vpb = vpb0./1000; % brine velocity in km/s
K_brine = vpb.*vpb.*rho_brine; % GPa