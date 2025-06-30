% Exact Reflection and Transmission Equations for VTI media
% see Graebner, 1992 and Rueger, 2002 and Luo, Ba and Carcione, 2020 papers
% in this version the inputs are velocities and Thomsen parameters
% coded by Ufuk Durmus on 2/2025

function [RP, RPS, thetap] = Graebner_FUN_vel(vp1,vs1,rho1,epsilon1,delta1,vp2,vs2,rho2,epsilon2,delta2)

% Proper Input Conversions
C33_1 = vp1^2*rho1;
C55_1 = vs1^2*rho1;
C11_1 = (2*epsilon1+1)*C33_1;
C13_1 = sqrt(2*delta1*C33_1*(C33_1-C55_1)+(C33_1-C55_1)^2)-C55_1;

C33_2 = vp2^2*rho2;
C55_2 = vs2^2*rho2;
C11_2 = (2*epsilon2+1)*C33_2;
C13_2 = sqrt(2*delta2*C33_2*(C33_2-C55_2)+(C33_2-C55_2)^2)-C55_2;

% the method is based on ray parameter not angle of incidence
%	--->>>  specify maximum ray parameter
%	      ------------------------------------
	pmax = 1  / vs1;   % reflection coefficients will be complex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	--->>>    linearly interpolate to get ray parameters
%		----------------------------------------------
	p = linspace( 0.0, pmax, 100)';

for i=1:length(p)
% for upper layer
% K1_1 = (rho1/C33_1)+(rho1/C55_1)-(C11_1/C55_1+C55_1/C33_1-(((C11_1+C55_1)^2)/(C33_1*C55_1)))*(p^2); % different from Rueger, 2002 C11 instead of C13, double CHECK!
K1_1 = (rho1/C33_1)+(rho1/C55_1)-(C11_1/C55_1+C55_1/C33_1-(((C13_1+C55_1)^2)/(C33_1*C55_1)))*(p(i)^2); % different from Rueger, 2002 C11 instead of C13, double CHECK!
K2_1 = (C11_1/C33_1)*(p(i)^2)-(rho1/C33_1);
K3_1 = (p(i)^2)-(rho1/C55_1);

sp1 = (1/sqrt(2))*sqrt(K1_1-sqrt(K1_1^2-4*K2_1*K3_1));
ss1 = (1/sqrt(2))*sqrt(K1_1+sqrt(K1_1^2-4*K2_1*K3_1));

lp1 = sqrt((C33_1*sp1^2+C55_1*p(i)^2-rho1)/...
    (C55_1*sp1^2+C11_1*p(i)^2-rho1+C33_1*sp1^2+C55_1*p(i)^2-rho1));

np1 = sqrt((C55_1*sp1^2+C11_1*p(i)^2-rho1)/...
    (C55_1*sp1^2+C11_1*p(i)^2-rho1+C33_1*sp1^2+C55_1*p(i)^2-rho1));

ls1 = sqrt((C55_1*ss1^2+C11_1*p(i)^2-rho1)/...
    (C55_1*ss1^2+C11_1*p(i)^2-rho1+C33_1*ss1^2+C55_1*p(i)^2-rho1));

ns1 = sqrt((C33_1*ss1^2+C55_1*p(i)^2-rho1)/...
    (C55_1*ss1^2+C11_1*p(i)^2-rho1+C33_1*ss1^2+C55_1*p(i)^2-rho1));

% for lower layer
% K1_2 = (rho2/C33_2)+(rho2/C55_2)-(C11_2/C55_2+C55_2/C33_2-(((C11_2+C55_2)^2)/(C33_2*C55_2)))*(p^2); % different from Rueger, 2002 C11 instead of C13, double CHECK!
K1_2 = (rho2/C33_2)+(rho2/C55_2)-(C11_2/C55_2+C55_2/C33_2-(((C13_2+C55_2)^2)/(C33_2*C55_2)))*(p(i)^2); % different from Rueger, 2002 C11 instead of C13, double CHECK!
K2_2 = (C11_2/C33_2)*(p(i)^2)-(rho2/C33_2);
K3_2 = (p(i)^2)-(rho2/C55_2);

sp2 = (1/sqrt(2))*sqrt(K1_2-sqrt(K1_2^2-4*K2_2*K3_2));
ss2 = (1/sqrt(2))*sqrt(K1_2+sqrt(K1_2^2-4*K2_2*K3_2));

lp2 = sqrt((C33_2*sp2^2+C55_2*p(i)^2-rho2)/...
    (C55_2*sp2^2+C11_2*p(i)^2-rho2+C33_2*sp2^2+C55_2*p(i)^2-rho2));

np2 = sqrt((C55_2*sp2^2+C11_2*p(i)^2-rho2)/...
    (C55_2*sp2^2+C11_2*p(i)^2-rho2+C33_2*sp2^2+C55_2*p(i)^2-rho2));

ls2 = sqrt((C55_2*ss2^2+C11_2*p(i)^2-rho2)/...
    (C55_2*ss2^2+C11_2*p(i)^2-rho2+C33_2*ss2^2+C55_2*p(i)^2-rho2));

ns2 = sqrt((C33_2*ss2^2+C55_2*p(i)^2-rho2)/...
    (C55_2*ss2^2+C11_2*p(i)^2-rho2+C33_2*ss2^2+C55_2*p(i)^2-rho2));

% S Matrix Variables
a1 = C55_1*(sp1*lp1+p(i)*np1);
b1 = C55_1*(ss1*ns1-p(i)*ls1);
a2 = C55_2*(sp2*lp2+p(i)*np2);
b2 = C55_2*(ss2*ns2-p(i)*ls2);
d1 = p(i)*lp1*C13_1+sp1*np1*C33_1;
e1 = p(i)*ns1*C13_1-ss1*ls1*C33_1;
d2_ = -p(i)*lp2*C13_2-sp2*np2*C33_2;
e2_ = -p(i)*ns2*C13_2+ss2*ls2*C33_2;

% S Matrix Elements
S(1,1) = lp1;
S(1,2) = ns1;
S(1,3) = -lp2;
S(1,4) = -ns2;
S(2,1) = np1;
S(2,2) = -ls1;
S(2,3) = np2;
S(2,4) = -ls2;
S(3,1) = a1;
S(3,2) = b1;
S(3,3) = a2;
S(3,4) = b2;
S(4,1) = d1;
S(4,2) = e1;
S(4,3) = d2_;
S(4,4) = e2_;

% B Matrix
B = [-lp1, np1, C55_1*(sp1*lp1+p(i)*np1), -p(i)*lp1*C13_1-sp1*np1*C33_1];
B = transpose(B);

% Reflection and Transmission Coefficeints RC = S^-1*B
RC = inv(S)*B;
RP(i) = RC(1,1);
RPS(i) = RC(2,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	--->>>   reflection coefficients may be complex depending upon
%		  the maximum ray parameter
%	--->>>   plot the real part: ie  real(pp),  real(ss) ... etc
%
%	--->>> to display as function of incidence angle (rather than ray parameter)
%		thetap = asin( p * alpha1 ) * 180 / pi ;
%		thetas = asin( p *  beta1 ) * 180 / pi ;
%
		thetap = asin( p * vp1 ) * 180 / pi ;
end