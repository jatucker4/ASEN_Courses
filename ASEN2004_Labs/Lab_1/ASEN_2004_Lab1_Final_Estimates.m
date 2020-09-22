%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear;
close all;
clc;
%% Extrapolating maximum height
[T, a, P, rho] = atmoscoesa(1624, 'None');
a = 0.0870;
alpha_L_0 = 0;
s_wing = 47564.0309263 / (1000^2) * 2; %[m^2]
t_crash = 5;
height = 6.4;
d = 42.6;
V_trim = sqrt(height^2 + d^2)/t_crash;
W = (78/1000)*9.81;
C_L = W/(0.5*rho*s_wing*V_trim^2);
alpha = (C_L/a)+alpha_L_0;
C_D = 0.0505;
L_D = C_L/C_D;