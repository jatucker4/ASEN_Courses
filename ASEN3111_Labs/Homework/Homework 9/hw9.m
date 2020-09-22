%% Code for 9.10
clc;clear all;close all;
M1 = 2; 
theta = 23.28*(pi/180);
p1 = .7;
T1 = 630;

v1 = pm_eqn_v(M1);
v1_deg = v1*(180/pi);

v2 = theta + v1;
M2 = pm_eqn_m(v2);

T2 = T1*((1+(.4/2)*M1^2)/(1+(.4/2)*M2^2));
p2 = p1*((1+(.4/2)*M1^2)/(1+(.4/2)*M2^2))^(1.4/.4);

rho2 = p2*2116.22/(1716*T2);

p02 = p2*(1+(.4/2)*M2^2)^(1.4/.4);
T02 = T2*(1+(.4/2)*M2^2);

mu1 = asind(1/M1);
mu2 = asind(1/M2);
aorma = mu2 - theta*(180/pi);

%% Begin code block for 9.13
clc; clear all; close all;
M1 = 2.6;
alpha = 5;
[cl_a,cd_a] = get_cl_cd_aoa(M1,alpha);

alpha = 15;
[cl_b,cd_b] = get_cl_cd_aoa(M1,alpha);

alpha = 30;
[cl_c,cd_c] = get_cl_cd_aoa(M1,alpha);



%% Begin code block for 9.14
clc;clear all;close all;
M1 = 3; 
theta2 = 5*(pi/180);
theta3 = 20*(pi/180);

v1 = pm_eqn_v(M1);
v1_deg = v1*(180/pi);

v2 = theta2 + v1;
v2_deg = v2*(180/pi);
M2 = -pm_eqn_m(v2);

p2op1 = pop(M1,M2);

v3 = v2 + theta3;
v3_deg = v3*(180/pi);
M3 = -pm_eqn_m(v3);

p3op1 = pop(M1,M3);

beta1 = get_beta(25*pi/180,M1);
Mn1 = M1*sin(beta1);
p4op1 = get_p2op1(Mn1);

Mn4 = get_Mn2(Mn1);
M4 = Mn4/sin(beta1-(25*pi/180));

v4 = pm_eqn_v(M4);
v4_deg = v4*180/pi;

v5 = v4 + theta3;
M5 = pm_eqn_m(v5);
p5op4 = pop(M4,M5);
p5op1 = p5op4*p4op1;

l_c = 1/(2*cosd(10));
prop = 2/(1.4*M1^2);
cl = prop*l_c*( (p4op1-p3op1)*cosd(25) + (p5op1-p2op1)*cosd(5));

cd = prop*l_c*( (p4op1-p3op1)*sind(25) + (p5op1-p2op1)*sind(5));

%% Functions used in this assignment
% <include> get_beta.m get_cl_cd_aoa.m get_Mn2.m get_p2op1.m pm_eqn_m.m
% pm_eqn_v.m pop.m </include>