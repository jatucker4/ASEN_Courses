%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is for the fourth ASEN3200 Homework
% 
% Created by: Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
%% Create constants
r_1 = [-1.054*10^8, 1.579*10^8, -1.520*10^5];
r_2 = [-1.461*10^8, 1.081*10^8, -2.265*10^5];
r_3 = [-1.652*10^8, 4.254*10^7, -2.673*10^5];
mu = 132712000000.000;

%% Use Gibb's function to calculate velocity
[r_out, v_out] = Gibbs_method(r_1,r_2,r_3,mu);

%% Use conversion function to get orbital elements
oevec = r_v_to_elements(r_2,v_out,mu);

%% Calcuate the time to cross the ecliptic
theta = oevec(7)*(pi/180);
e = oevec(3);
a = oevec(2);
ohm = 2*pi - oevec(5)*(pi/180);

E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));
E = (2*pi) + E;
Me = E - e*sin(E);
T = sqrt((a^3)/mu)*(2*pi);
t = (Me/(2*pi))*T;

E_1 = 2*atan(sqrt((1-e)/(1+e))*tan(ohm/2));
Me_1 = E_1 - e*sin(E_1);
T_1 = sqrt((a^3)/mu)*(2*pi);
t_1 = (Me_1/(2*pi))*T_1;

t_total = T - t + t_1;