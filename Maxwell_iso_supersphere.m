function [K_eff,mu_eff]= Maxwell_iso_supersphere(K0, mu0, s, por)
% based on the paper Chen, Sevostianov, Giraud and Grgic, 2015 and see also Sevostianov and Giraud, 2013
% coded by Ufuk Durmus on 11/2024
% INPUTS
% K0 =  Bulk modulus of the background matrix
% mu0 =  Shear modulus of the background matrix
% s = non-spherity parameter <0.5 concave, >0.5 convex, =1 sphere shape
% por =  volume fraction of inclusion
% OUTPUTS; effective isotropic bulk (K) and shear (mu) moduli
%% Input here

lambda0 = K0-(2/3)*mu0; % Lame parameter of the background matrix
v0 = lambda0/(2*(lambda0+mu0));  % poisson's ratio of the background matrix

%% STEPS
% phi parameters for bulk and shear moduli calculation
phi_K = (2/3)*(1-2*v0)/(1-v0);
phi_mu = (1/15)*(7-5*v0)/(1-v0);

% see Chen, Sevostianov, Giraud and Grgic, 2015 defined differently
H_K = (1-v0)/(1-2*v0);
H_mu = 10*(1-v0)/(7-5*v0);

% x microstructural parameter depends on s non-spherity parameter
% s is defined as a supersphere in x^2s+y^2s+z^2s=1 
% gamma is defined as integral, see the papers
fun1 = @(t) exp(-t).*(t.^((3/(2.*s))-1));
gamma1 = integral(fun1,0,Inf);
fun2 = @(t) exp(-t).*(t.^((1/(2.*s))-1));
gamma2 = integral(fun2,0,Inf);

% see Chen, Sevostianov, Giraud and Grgic, 2015 defined differently
x = (3*pi/4*(5*s-1)*(s^2)*gamma1)/(gamma2^3);

%% LAST STEP: effective elastic properties estimation in an isotropic medium
K_ratio = (1-por.*x.*H_K.*phi_K)./(1+por.*x.*H_K.*(1-phi_K));
mu_ratio = (1-por.*x.*H_mu.*phi_mu)./(1+por.*x.*H_mu.*(1-phi_mu));

K_eff = K_ratio.*K0;
mu_eff = mu_ratio.*mu0;

end