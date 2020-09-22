function [cl,cd] = get_cl_cd_aoa(M1,alpha)
%GET_CL_CD_AOA Summary of this function goes here
%   Detailed explanation goes here
alpha = alpha*pi/180;
v1 = pm_eqn_v(M1);
v1_deg = v1*180/pi;

v2 = v1 + alpha;
M2 = pm_eqn_m(v2);
p2op1 = pop(M1,M2);

beta1 = get_beta(alpha,M1);
Mn1 = M1*sin(beta1);
p3op1 = 1 + (2.8/2.4)*(Mn1^2 - 1);

cl = (2/(1.4*M1^2))*(p3op1 - p2op1)*cos(alpha);
cd = (2/(1.4*M1^2))*(p3op1 - p2op1)*sin(alpha);
end

