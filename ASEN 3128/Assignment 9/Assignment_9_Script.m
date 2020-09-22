clear;close all; clc;

rho = 0.3045;
b = 59.64;
S = 511;
u_0 = 235.9;

C_y_B = -0.8771;
C_y_p = 0;
C_y_r = 0;

C_l_B = -0.2797;
C_l_p = -0.3295;
C_l_r = 0.304;

C_n_B = 0.1964;
C_n_p = -0.04073;
C_n_r = -0.2737;

Y_v = 0.5*rho*u_0*S*C_y_B;
Y_p = 0.25*rho*u_0*b*S*C_y_p;
Y_r = 0.25*rho*u_0*b*S*C_y_r;

L_v = 0.5*rho*u_0*b*S*C_l_B;
L_p = 0.25*rho*u_0*(b^2)*S*C_l_p;
L_r = 0.25*rho*u_0*(b^2)*S*C_l_r;

N_v = 0.5*rho*u_0*b*S*C_n_B;
N_p = 0.25*rho*u_0*(b^2)*S*C_n_p;
N_r = 0.25*rho*u_0*(b^2)*S*C_n_r;