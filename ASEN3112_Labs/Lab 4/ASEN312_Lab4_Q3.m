%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3112 - Lab 4
% 
% Created By: Johnathan Tucker
%
% Collaborators: N/A
%
% The purpose of the script is to contain all of the constants and basic
% computation that will be used to answer the third question in Lab 4
%
% Created Date: 4/23/2020
%
% Change Log: 
%           - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
tic
%Global formatting commands to imporve graphing looks:
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex'); 
%% Create necessary constants and equations
E = 10e6; %[psi]
YS = 35e3; %[psi]
% I_1 = (.0625^3)*.25/12; %[in^4]
I_1 = 1.63e-4;
% I_2 = (.125^3)/12; %[in^4]
I_2 = 3e-4;
A_1 = .25^2 - (.25-2*.0625)^2; %[in^2]
A_2 = .125*1; %[in^2]
% Length vector
L = linspace(2,12); % [in]
%% Calculate the buckling load for each case
% First calculate pinned-pinned buckling for each
Pcr_pp_1 = ((pi.^2).*E.*I_1)./(L.^2); % [lbf]
Pcr_pp_2 = ((pi.^2).*E.*I_2)./(L.^2); % [lbf]

% Next calculate the fixed-fixed buckling for each
Pcr_ff_1 = ((pi.^2).*E.*I_1)./((.5.*L).^2); % [lbf]
Pcr_ff_2 = ((pi.^2).*E.*I_2)./((.5.*L).^2); % [lbf]

% Now calculate the yield Load for each beam
YL_1 = YS.*A_1; % [lbf]
YL_2 = YS.*A_2; % [lbf]
%% Create Plots
figure
hold on
plot(L, Pcr_pp_1)
plot(L, Pcr_ff_1)
plot(L, ones(1,length(L)).*YL_1)
hold off
legend("Pinned-Pinned ","Fixed-Fixed","Yield Load")
xlabel("Length [in]")
ylabel("Buckling Load [lbf]")
title("Buckling Load vs Beam Length for Hollow Square Rod")

figure
hold on
plot(L, Pcr_pp_2)
plot(L, Pcr_ff_2)
plot(L, ones(1,length(L)).*YL_2)
hold off
legend("Pinned-Pinned ","Fixed-Fixed","Yield Load")
xlabel("Length [in]")
ylabel("Buckling Load [lbf]")
title("Buckling Load vs Beam Length for Solid Rectangular Rod")


