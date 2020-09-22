%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is for problem 18-19 in homework 7
%
% Created by: Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
%% Create variables
L_c = 0.01;
rho = 8000;
c_p = 570;
T_0 = 291;
T_inf = 1223;
h = 150;
d = 3;
v_range = linspace(0.005,0.060);
t = d./v_range;
b = (h/(rho*c_p*L_c));
%% Calculate the temperatures and plot
T_vec = T_inf + exp(-b.*t).*(T_0-T_inf);
T_vec = T_vec - 273;

plot(v_range,T_vec)
xlabel("$Velocity\:Range\:[m/s]$",'Interpreter','latex','FontSize',26)
ylabel("$Plate\:Temperature\:[C]$",'Interpreter','latex','FontSize',26)
title("$Plate\:Temperature\:vs\:Velocity$",...
    'Interpreter','latex','FontSize',26)