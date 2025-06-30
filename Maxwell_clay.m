function [C11,C33,C13,C55,C66]= Maxwell_clay (K0,G0,K1,G1,K2,G2,K3,G3,v1,v2,v3,gamma1,gamma2,gamma3)

% Coded by Ufuk Durmus
% Last Updated 3/16/2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs here

G = G0;
mu = G0;

lambda0 = K0-(2/3)*G0; %Lame Parameter of the matrix

v0 = lambda0/(2*(lambda0+G0));  % poisson's ratio of the matrix

lambda1 = K1-(2/3)*G1; %Lame Parameter of inclusion 1

lambda2 = K2-(2/3)*G2; %Lame Parameter of inclusion 2

lambda3 = K3-(2/3)*G3; %Lame Parameter of inclusion 3

% Define elastic stiffness tensor of the matrix

Cij_0=zeros(6,6);
for i=1:6
    if i<4
        Cij_0(i,i)=lambda0+2*G0;
    else
        Cij_0(i,i)=G0;
    end
    Cij_0(1,2)=lambda0;
    Cij_0(1,3)=lambda0;
    Cij_0(2,1)=lambda0;
    Cij_0(3,1)=lambda0;
    Cij_0(3,2)=lambda0;
    Cij_0(2,3)=lambda0;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion Using Tensorial Basis

C0_1 = (Cij_0(1,1)+Cij_0(1,2))/2;
C0_2 = 2*Cij_0(6,6);
C0_3 = Cij_0(1,3);
C0_4 = Cij_0(3,1);
C0_5 = 4*Cij_0(5,5);
C0_6 = Cij_0(3,3);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Inclusion 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gamma1 < 1 % oblate shape
g_gamma_1 = (1/(gamma1*sqrt(1-gamma1^2)))*atan(sqrt(1-gamma1^2)/gamma1);
elseif gamma1 > 1 % prolate shape
    g_gamma_1 = (1/(2*gamma1*sqrt(gamma1^2-1)))*(log((gamma1+sqrt(gamma1^2-1))...
        /(gamma1-sqrt(gamma1^2-1))));
else
    error('based on the paper, gamma should be either greater or lower than 1, modify it!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1/(2*(1-v0));
f0_1 = (gamma1^2*(1-g_gamma_1))/(2*(gamma1^2-1));
f1_1 = ((k*gamma1^2)/(4*(gamma1^2-1)^2))*((2*gamma1^2+1)*g_gamma_1-3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Inclusion 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gamma2 < 1 % oblate shape
g_gamma_2 = (1/(gamma2*sqrt(1-gamma2^2)))*atan(sqrt(1-gamma2^2)/gamma2);
elseif gamma2 > 1 % prolate shape
    g_gamma_2 = (1/(2*gamma2*sqrt(gamma2^2-1)))*(log((gamma2+sqrt(gamma2^2-1))...
        /(gamma2-sqrt(gamma2^2-1))));
else
    error('based on the paper, gamma should be either greater or lower than 1, modify it!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1/(2*(1-v0));
f0_2 = ((gamma2^2)*(1-g_gamma_2))/(2*(gamma2^2-1));
f1_2 = ((k*(gamma2^2))/(4*(gamma2^2-1)^2))*((2*gamma2^2+1)*g_gamma_2-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Inclusion 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gamma3 < 1 % oblate shape
g_gamma_3 = (1/(gamma3*sqrt(1-gamma3^2)))*atan(sqrt(1-gamma3^2)/gamma3);
elseif gamma3 > 1 % prolate shape
    g_gamma_3 = (1/(2*gamma3*sqrt(gamma3^2-1)))*(log((gamma3+sqrt(gamma3^2-1))...
        /(gamma3-sqrt(gamma3^2-1))));
else
    error('based on the paper, gamma should be either greater or lower than 1, modify it!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1/(2*(1-v0));
f0_3 = ((gamma3^2)*(1-g_gamma_3))/(2*(gamma3^2-1));
f1_3 = ((k*(gamma3^2))/(4*(gamma3^2-1)^2))*((2*gamma3^2+1)*g_gamma_3-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The stiffness contribution tensor N for a spheroidal inhomogeneity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inclusion 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1_1 = (1/(2*mu))*((1-k)*f0_1+k*f1_1);
p2_1 = (1/(2*mu))*((2-k)*f0_1+k*f1_1);
p3_1 = - (k/mu)*f1_1;
p4_1 = - (k/mu)*f1_1;
p5_1 = (1/mu)*(1-f0_1-4*k*f1_1);
p6_1 = (1/mu)*((1-k)*(1-2*f0_1)+(2*k*f1_1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_G_n_1 = G1 - G0;
delta_lambda_n_1 = lambda1 - lambda0;
delta1_n_1 = ((1+(delta_lambda_n_1+2*delta_G_n_1)*p6_1+4*(delta_lambda_n_1+delta_G_n_1)*p1_1+4*delta_lambda_n_1*p3_1)...
    /(2*delta_G_n_1*(3*delta_lambda_n_1+2*delta_G_n_1)))+(2*p1_1*p6_1-2*p3_1^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1_1 = (1/(2*delta1_n_1))*((delta_lambda_n_1+delta_G_n_1)/(delta_G_n_1*(3*delta_lambda_n_1+2*delta_G_n_1))...
    +p6_1);
n2_1 = (2*delta_G_n_1)/(1+2*p2_1*delta_G_n_1);
n3_1 = (-1/delta1_n_1)*((-delta_lambda_n_1)/(2*delta_G_n_1*(3*delta_lambda_n_1+2*delta_G_n_1))+p3_1);
n4_1 = (-1/delta1_n_1)*((-delta_lambda_n_1)/(2*delta_G_n_1*(3*delta_lambda_n_1+2*delta_G_n_1))+p3_1);
n5_1 = (4*delta_G_n_1)/(1+delta_G_n_1*p5_1);
n6_1 = (1/delta1_n_1)*((delta_lambda_n_1+2*delta_G_n_1)/(2*delta_G_n_1*(3*delta_lambda_n_1+2*delta_G_n_1))+2*p1_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inclusion 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1_2 = (1/(2*mu))*((1-k)*f0_2+k*f1_2);
p2_2 = (1/(2*mu))*((2-k)*f0_2+k*f1_2);
p3_2 = - (k/mu)*f1_2;
p4_2 = - (k/mu)*f1_2;
p5_2 = (1/mu)*(1-f0_2-4*k*f1_2);
p6_2 = (1/mu)*((1-k)*(1-2*f0_2)+(2*k*f1_2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_G_n_2 = G2 - G0;
delta_lambda_n_2 = lambda2 - lambda0;
delta1_n_2 = ((1+(delta_lambda_n_2+2*delta_G_n_2)*p6_2+4*(delta_lambda_n_2+delta_G_n_2)*p1_2+4*delta_lambda_n_2*p3_2)...
    /(2*delta_G_n_2*(3*delta_lambda_n_2+2*delta_G_n_2)))+(2*p1_2*p6_2-2*(p3_2^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1_2 = (1/(2*delta1_n_2))*((delta_lambda_n_2+delta_G_n_2)/(delta_G_n_2*(3*delta_lambda_n_2+2*delta_G_n_2))...
    +p6_2);
n2_2 = (2*delta_G_n_2)/(1+2*p2_2*delta_G_n_2);
n3_2 = (-1/delta1_n_2)*((-delta_lambda_n_2)/(2*delta_G_n_2*(3*delta_lambda_n_2+2*delta_G_n_2))+p3_2);
n4_2 = (-1/delta1_n_2)*((-delta_lambda_n_2)/(2*delta_G_n_2*(3*delta_lambda_n_2+2*delta_G_n_2))+p3_2);
n5_2 = (4*delta_G_n_2)/(1+delta_G_n_2*p5_2);
n6_2 = (1/delta1_n_2)*((delta_lambda_n_2+2*delta_G_n_2)/(2*delta_G_n_2*(3*delta_lambda_n_2+2*delta_G_n_2))+2*p1_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inclusion 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1_3 = (1/(2*mu))*((1-k)*f0_3+k*f1_3);
p2_3 = (1/(2*mu))*((2-k)*f0_3+k*f1_3);
p3_3 = - (k/mu)*f1_3;
p4_3 = - (k/mu)*f1_3;
p5_3 = (1/mu)*(1-f0_3-4*k*f1_3);
p6_3 = (1/mu)*((1-k)*(1-2*f0_3)+(2*k*f1_3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_G_n_3 = G3 - G0;
delta_lambda_n_3 = lambda3 - lambda0;
delta1_n_3 = ((1+(delta_lambda_n_3+2*delta_G_n_3)*p6_3+4*(delta_lambda_n_3+delta_G_n_3)*p1_3+4*delta_lambda_n_3*p3_3)...
    /(2*delta_G_n_3*(3*delta_lambda_n_3+2*delta_G_n_3)))+(2*p1_3*p6_3-2*(p3_3^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1_3 = (1/(2*delta1_n_3))*((delta_lambda_n_3+delta_G_n_3)/(delta_G_n_3*(3*delta_lambda_n_3+2*delta_G_n_3))...
    +p6_3);
n2_3 = (2*delta_G_n_3)/(1+2*p2_3*delta_G_n_3);
n3_3 = (-1/delta1_n_3)*((-delta_lambda_n_3)/(2*delta_G_n_3*(3*delta_lambda_n_3+2*delta_G_n_3))+p3_3);
n4_3 = (-1/delta1_n_3)*((-delta_lambda_n_3)/(2*delta_G_n_3*(3*delta_lambda_n_3+2*delta_G_n_3))+p3_3);
n5_3 = (4*delta_G_n_3)/(1+delta_G_n_3*p5_3);
n6_3 = (1/delta1_n_3)*((delta_lambda_n_3+2*delta_G_n_3)/(2*delta_G_n_3*(3*delta_lambda_n_3+2*delta_G_n_3))+2*p1_3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating Q and P tensors for the shape of omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inclusion 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q1_1 = G*(4*k-1-2*(3*k-1)*f0_1-2*k*f1_1);
q2_1 = 2*G*(1-(2-k)*f0_1-k*f1_1);
q3_1= 2*G *((2*k-1)*f0_1+2*k*f1_1);
q4_1 = 2*G*((2*k-1)*f0_1+2*k*f1_1);
q5_1 = 4*G*(f0_1+4*k*f1_1);
q6_1 = 8*G*(k*f0_1-k*f1_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inclusion 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q1_2 = G*(4*k-1-2*(3*k-1)*f0_2-2*k*f1_2);
q2_2 = 2*G*(1-(2-k)*f0_2-k*f1_2);
q3_2 = 2*G *((2*k-1)*f0_2+2*k*f1_2);
q4_2 = 2*G*((2*k-1)*f0_2+2*k*f1_2);
q5_2 = 4*G*(f0_2+4*k*f1_2);
q6_2 = 8*G*(k*f0_2-k*f1_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inclusion 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q1_3 = G*(4*k-1-2*(3*k-1)*f0_3-2*k*f1_3);
q2_3 = 2*G*(1-(2-k)*f0_3-k*f1_3);
q3_3 = 2*G *((2*k-1)*f0_3+2*k*f1_3);
q4_3 = 2*G*((2*k-1)*f0_3+2*k*f1_3);
q5_3 = 4*G*(f0_3+4*k*f1_3);
q6_3 = 8*G*(k*f0_3-k*f1_3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q_1(3,3) = (q6_1);
Q_1(1,3) = (q3_1);
Q_1(5,5) = (q5_1./4);
Q_1(6,6) = (q2_1./2);
Q_1(1,1) = (q1_1+Q_1(6,6));

Q_2(3,3) = (q6_2);
Q_2(1,3) = (q3_2);
Q_2(5,5) = (q5_2./4);
Q_2(6,6) = (q2_2./2);
Q_2(1,1) = (q1_2+Q_2(6,6));

Q_3(3,3) = (q6_3);
Q_3(1,3) = (q3_3);
Q_3(5,5) = (q5_3./4);
Q_3(6,6) = (q2_3./2);
Q_3(1,1) = (q1_3+Q_3(6,6));

% Q_4(3,3) = (q6_4);
% Q_4(1,3) = (q3_4);
% Q_4(5,5) = (q5_4./4);
% Q_4(6,6) = (q2_4./2);
% Q_4(1,1) = (q1_4+Q_4(6,6));

omega_q = (v1*Q_1(3,3)+v2*Q_2(3,3)+v3*Q_3(3,3))...
/(v1*Q_1(1,1)+v2*Q_2(1,1)+v3*Q_3(1,1));


% P_1(3,3) = (p6_1);
% P_1(1,3) = (p3_1);
% P_1(5,5) = (p5_1./4);
% P_1(6,6) = (p2_1./2);
% P_1(1,1) = (p1_1+P_1(6,6));
% 
% P_2(3,3) = (p6_2);
% P_2(1,3) = (p3_2);
% P_2(5,5) = (p5_2./4);
% P_2(6,6) = (p2_2./2);
% P_2(1,1) = (p1_2+P_2(6,6));
% 
% omega_q = (v1*P_1(1,1)+v2*P_2(1,1))/(v1*P_1(3,3)+v2*P_2(3,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if omega_q < 1 % oblate shape
g_omega_q = (1/(omega_q*sqrt(1-omega_q^2)))*atan(sqrt(1-omega_q^2)/omega_q);
elseif omega_q > 1 % prolate shape
    g_omega_q = (1/(2*omega_q*sqrt(omega_q^2-1)))*(log((omega_q+sqrt(omega_q^2-1))...
        /(omega_q-sqrt(omega_q^2-1))));
else
    error('based on the paper, gamma should be either greater or lower than 1, modify it!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1/(2*(1-v0));
f0_omega_q = (omega_q.^2*(1-g_omega_q))/(2*(omega_q.^2-1));
f1_omega_q = ((k*omega_q^2)/(4*(omega_q^2-1)^2))*((2*omega_q^2+1)*g_omega_q-3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q1 = G*(4*k-1-2*(3*k-1)*f0_omega_q-2*k*f1_omega_q);
Q2 = 2*G*(1-(2-k)*f0_omega_q-k*f1_omega_q);
Q3 = 2*G *((2*k-1)*f0_omega_q+2*k*f1_omega_q);
Q4 = 2*G*((2*k-1)*f0_omega_q+2*k*f1_omega_q);
Q5 = 4*G*(f0_omega_q+4*k*f1_omega_q);
Q6 = 8*G*(k*f0_omega_q-k*f1_omega_q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimating P tensor

% P = S_0*(J-Q*S_0)
% Let's say X = Q*S_0
% and let's say M = (J-X)
% P = S_0*M

% Converting stiffness tensor to compliance tensor

S_1 = (1-v0)/(4*G*(1+v0));
S_2 = 1/(2*G);
S_3 = -(v0)/(2*G*(1+v0));
S_4 = -(v0)/(2*G*(1+v0));
S_5 = 1/G;
S_6 = 1/(2*G*(1+v0));

J1 = 1/2;
J2 = 1;
J3 = 0;
J4 = 0;
J5 = 2;
J6 = 1;

X1 = (2*Q1*S_1+Q3*S_4);
X2 = (Q2*S_2);
X3 = (2*Q1*S_3+Q3*S_6);
X4 = (2*Q4*S_1+Q6*S_4);
X5 = (1/2)*(Q5*S_5);
X6 = (Q6*S_6+2*Q4*S_3);

M1 = J1-X1;
M2 = J2-X2;
M3 = J3-X3;
M4 = J4-X4;
M5 = J5-X5;
M6 = J6-X6;

P1 = (2*S_1*M1+S_3*M4);
P2 = (S_2*M2);
P3 = (2*S_1*M3+S_3*M6);
P4 = (2*S_4*M1+S_6*M4);
P5 = (1/2)*(S_5*M5);
P6 = (S_6*M6+2*S_4*M3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimating Effective Stiffness Tensor of the Heterogeneous Medium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C_eff = C_0 + inv(inv(v1.*N_1+v2.*N_2)-P)
% let's say N = inv(v1.*N_1+v2.*N_2)
% and let's say T = inv(N-P)
% and C_eff = C_0 + T

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N1 = v1*n1_1+v2*n1_2+v3*n1_3;
N2 = v1*n2_1+v2*n2_2+v3*n2_3;
N3 = v1*n3_1+v2*n3_2+v3*n3_3;
N4 = v1*n4_1+v2*n4_2+v3*n4_3;
N5 = v1*n5_1+v2*n5_2+v3*n5_3;
N6 = v1*n6_1+v2*n6_2+v3*n6_3;

% Inverse of N using tensorial basis
delta_n = 2*(N1*N6-N3*N4);

N1_i = (N6/(2*delta_n));
N2_i = (1/N2);
N3_i = -(N3/(delta_n));
N4_i = -(N4/(delta_n));
N5_i = (4/N5);
N6_i = ((2*N1)/(delta_n));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T1 = N1_i-P1;
T2 = N2_i-P2;
T3 = N3_i-P3;
T4 = N4_i-P4;
T5 = N5_i-P5;
T6 = N6_i-P6;

% Inverse of T using tensorial basis
delta_t = 2*(T1*T6-T3*T4);

T1_i = (T6/(2*delta_t));
T2_i = (1/T2);
T3_i = -(T3/(delta_t));
T4_i = -(T4/(delta_t));
T5_i = (4/T5);
T6_i = ((2*T1)/(delta_t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C1 = C0_1+T1_i;
C2 = C0_2+T2_i;
C3 = C0_3+T3_i;
C4 = C0_4+T4_i;
C5 = C0_5+T5_i;
C6 = C0_6+T6_i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cij_eff=zeros(6,6);

Cij_eff(3,3) = (C6);
Cij_eff(1,3) = (C3);
Cij_eff(5,5) = (C5./4);
Cij_eff(6,6) = (C2./2);
Cij_eff(1,1) = (C1+Cij_eff(6,6));
Cij_eff(2,2) = (Cij_eff(1,1));
Cij_eff(2,2) = (Cij_eff(1,1));
Cij_eff(3,1) = (Cij_eff(1,3));
Cij_eff(3,2) = (Cij_eff(1,3));
Cij_eff(2,3) = (Cij_eff(1,3));
Cij_eff(4,4) = (Cij_eff(5,5));
Cij_eff(2,1) = (Cij_eff(1,1)-2*Cij_eff(6,6));
Cij_eff(1,2) = (Cij_eff(1,1)-2*Cij_eff(6,6));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Five Independent Stiffness Coefficient
C11 = Cij_eff(1,1);
C33 = Cij_eff(3,3);
C13 = Cij_eff(1,3);
C55 = Cij_eff(5,5);
C66 = Cij_eff(6,6);

end