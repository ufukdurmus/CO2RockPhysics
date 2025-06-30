function[R_vti_p,R_vti_ps,R_iso_ps,R_vti_sv,R_vti_sh]=rueger_vti_avo_fun(a1,b1,rho1,a2,b2,rho2,e1,g1,d1,e2,g2,d2,theta)
% Rueger - calculate the P, PS, SV and SH-wave reflectivity 
% in weakly anisotropic VTI media using Ruger's approximation
%
% a1, b1: P- and S-wave velocities of upper medium perpendicular to
% symmetry axis (vertical velocities for VTI)
% e1, g1, d1: epsilon gamma delta, Thomsen weakly anisotropic parameters of the upper medium
% a2, b2: P- and S-wave velocities perpendicular to symmetry axis in the lower medium
% e2,g2,d2: thomsen's aniso parameters of lower medium
% rho1, rho2: density of the upper and lower halfspace, respectively
% theta: the incidence angle
%
% R reflectivity as a function of incidence angle
%
% reference:
% Ruger, A., 1997, Geophysics

%Coded by Ufuk Durmus
%Last Updated on 3/16/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define variables

% %Test Rueger's Parameters
% 
% a1 = 3.30;
% b1 = 1.7;
% e1 = 0.133;
% g1 = 0.0;
% d1 = 0.12;
% rho1 = 2.35;
% 
% a2 = 4.20;
% b2 = 2.7;
% e2 = 0;
% g2 = 0;
% d2 = 0;
% rho2 = 2.49;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_z = (a2.*rho2-a1.*rho1);
avg_z = (a2.*rho2+a1.*rho1)./2;

delta_z_s = (b2.*rho2-b1.*rho1);
avg_z_s = (b2.*rho2+b1.*rho1)./2;

delta_alpha = (a2-a1);
avg_alpha = (a2+a1)./2;

delta_beta = (b2-b1);
avg_beta = (b2+b1)./2;

delta_rho = rho2-rho1;
avg_rho = (rho2+rho1)./2;

delta_G = ((b2.^2).*rho2-(b1.^2).*rho1);
avg_G = ((b2.^2).*rho2+(b1.^2).*rho1)./2;

delta_d = (d2-d1);
delta_e = (e2-e1);
delta_g = (g2-g1);
delta_dd = (d1-d2);
delta_ee = (e1-e2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:length(a1)
      
for i=1:length(theta)
        
R_vti_p(i,j) = (0.5*(delta_z(j)/avg_z(j)))+0.5*((delta_alpha(j)/avg_alpha(j))-((2*avg_beta(j))/avg_alpha(j))^2 ...
*(delta_G(j)/avg_G(j))+(delta_d(j)))*(sind(i))^2+0.5*(delta_alpha(j)/avg_alpha(j)+delta_e(j))*(sind(i))^2*(tand(i))^2;

R_iso_ps(i,j) = -0.5*(delta_rho(j)/avg_rho(j))*(sind(i)/cosd(i))-(avg_beta(j)/avg_alpha(j))...
    *((delta_rho(j)/avg_rho(j))+2*(delta_beta(j)/avg_beta(j)))*sind(i)*cosd(i)...
    +((avg_beta(j)/avg_alpha(j))^2)*(2*(delta_beta(j)/avg_beta(j))+(delta_rho(j)/avg_rho(j)))...
    *(((sind(i))^3)/cosd(i));

R_vti_ps(i,j) = R_iso_ps(i,j) + (((avg_alpha(j)^2/(2*(avg_alpha(j)^2-avg_beta(j)^2)*cosd(i))...
    -(avg_alpha(j)*avg_beta(j)*cosd(i))/(2*(avg_alpha(j)^2-avg_beta(j)^2)))*delta_d(j))*sind(i))...
    +(((avg_beta(j)*avg_alpha(j)*cosd(i))/(avg_alpha(j)^2-avg_beta(j)^2))*(delta_d(j)+(delta_ee(j)))*(sind(i))^3)...
    -((avg_alpha(j)^2/((avg_alpha(j)^2-avg_beta(j)^2)*cosd(i)))*(delta_d(j)+(delta_ee(j)))*(sind(i))^3)...
    +((avg_beta(j)^2/(2*(avg_alpha(j)^2-avg_beta(j)^2)*cosd(i)))*(delta_dd(j))*(sind(i))^3)...
    +((avg_beta(j)^2/((avg_alpha(j)^2-avg_beta(j)^2)*cosd(i)))*(delta_d(j)+(delta_ee(j)))*(sind(i))^5);

R_vti_sv(i,j) = - 0.5*(delta_z_s(j)/avg_z_s(j))+((7/2)*(delta_beta(j)/avg_beta(j))+2*(delta_rho(j)/avg_rho(j))...
    +0.5*((avg_alpha(j)/avg_beta(j))^2)*(delta_e(j)-delta_d(j)))*(sind(i))^2 ...
-0.5*(delta_beta(j)/avg_beta(j))*(sind(i))^2*(tand(i))^2;

R_vti_sh(i,j) = - 0.5*(delta_z_s(j)/avg_z_s(j))+0.5*(delta_beta(j)/avg_beta(j)+delta_g(j))*(tand(i))^2;

end

end



end