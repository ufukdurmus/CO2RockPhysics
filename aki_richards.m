function [RC, I, G] = aki_richards(a1, b1, rho1, a2, b2, rho2, theta)
delta_alpha = (a2-a1); %p-wave velocity difference
avg_alpha = (a2+a1)./2; %p-wave velocity average

delta_beta = (b2-b1); %s-wave velocity difference
avg_beta = (b2+b1)./2; %s-wave velocity average

delta_rho = rho2-rho1; %density difference
avg_rho = (rho2+rho1)./2; %density average

gamma = (avg_beta./avg_alpha).^2; %Vs./Vp parameter

% Reflection Coefficient vs Angle based on Aki-Richards 3-term
for i = 1:length(theta)
RC(i) = (1/2).*(1-4.*gamma.*(sind(theta(i))).^2).*(delta_rho./avg_rho)...
    + (1./(2.*(cosd(theta(i))).^2)).*(delta_alpha./avg_alpha)...
    - (4.*gamma.*(sind(theta(i))).^2).*(delta_beta./avg_beta);
    
end

% Intercept & Gradient
I = 0.5.*((delta_alpha./avg_alpha)+(delta_rho./avg_rho));
G = I - 2.*(0.5.*((delta_beta./avg_beta)+(delta_rho./avg_rho)));

end