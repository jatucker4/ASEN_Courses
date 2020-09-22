%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is the driver for the first ASEN 3112 Lab
%
% Created by: Johnathan Tucker
%
% Created on: 2/13/2020
%
% Change Log: 
%               -2/13/2020: Initial commit. 
%               -2/16/2020: Corrected least squares theory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
close all;
clear all;

%% Begin Analysis of CTW section

% Create relevant CTW constants
t = 1/16; % [in]
R_e = 3/8; % [in]
R_i = R_e - t; % [in]
R_mean = (R_e + R_i)/2; % [in]
G = 3.75e6; % [psi]

% Load in the data
CTW_data = readmatrix('CTW.csv');

% Parse the data into variables relevant for analysis
gamma = CTW_data(:,3).*(pi/180); % shear strain [rad]
phi = CTW_data(:,2).*(pi/180); % twist angle [rad]
T = CTW_data(:,4); % Applied torque [lbs-in]
L = 9; % [in]
% Find where the torque hits its max(so we'll fit to the way down)
[~,index] = max(T);

% Use CTW theory for shear strain
shear_strain_ctw = (phi-phi(1))*R_e/L; % [rad]

% Perform least squares for each
[least_squares_ext, S_ext] = polyfit(T(1:index,1),gamma(1:index,1),1);
[least_squares_ctw, S_ctw] = polyfit(T(1:index,1),shear_strain_ctw(1:index,1),1);
GJ_ext_estimate = -R_mean/least_squares_ext(1);
GJ_ctw_estimate = -R_e/least_squares_ctw(1);

% Obtain the estimated error for each estimate
[~,ext_error] = polyval(least_squares_ext,T(1:index,1),S_ext);
[~,ctw_error] = polyval(least_squares_ctw,T(1:index,1),S_ctw);

% Calculate theoretical GJ's
exact_GJ = G*(pi/2)*(R_e^4 - R_i^4);
ctw_GJ = G*((4*(pi*(R_mean^2))^2)/(2*pi*R_e/t));

% Calculate percent the various percent differences
perc_diff_1 = (abs(exact_GJ - GJ_ext_estimate)/GJ_ext_estimate) * 100;
fprintf('The percent difference between the estimated rigidity using the extensometer data and the rigidy using exact theory is %f\n',perc_diff_1);
perc_diff_2 = (abs(exact_GJ - GJ_ctw_estimate)/GJ_ctw_estimate) * 100;
fprintf('The percent difference between the estimated rigidity using the twist angle and the rigidy using exact theory is %f\n',perc_diff_2);
perc_diff_3 = (abs(ctw_GJ - GJ_ext_estimate)/GJ_ext_estimate) * 100;
fprintf('The percent difference between the estimated rigidity using the extensometer data and the rigidy using ctw theory is %f\n',perc_diff_3);
perc_diff_4 = (abs(ctw_GJ - GJ_ctw_estimate)/GJ_ctw_estimate) * 100;
fprintf('The percent difference between the estimated rigidity using the twist angle data and the rigidy using ctw theory is %f\n\n',perc_diff_4);

% Create the necessary plots
figure(1)
plot(T,gamma)
hold on
plot(T,shear_strain_ctw)
legend("$Extensometer\:Data$","$Twist\:Angle\:Based$",'Interpreter','latex')
xlabel("$Torque\:[lbf-in]$",'Interpreter','latex','FontSize',26)
ylabel("$Shear\:Strain[\frac{in}{in}]$",'Interpreter','latex','FontSize',26)
title("$CTW\:Section\:Shear\:Strain\:vs\:Torque$",'Interpreter','latex','FontSize',26)
%% Begin analysis for the open thin walled section
% Create relevant OTW constants
t = 1/16; % [in]
R_e = 3/8; % [in]
R_i = R_e - t; % [in]
R_mean = (R_e + R_i)/2; % [in]
G = 3.75e6; % [psi]

% Load in the data
OTW_data = readmatrix('OTW.csv');

% Parse the data into variables relevant for analysis
gamma = (OTW_data(:,3)).*(pi/180); % shear strain [rad]
phi = (OTW_data(:,2)).*(pi/180); % twist angle [rad]
T = (OTW_data(:,4)); % Applied torque [lbs-in]
L = 9; % [in]
% Find where the torque hits its max(so we'll fit to the way down)
[~,index] = max(T);

% Use OTW theory for shear strain
shear_strain_otw = (phi-phi(1))*t/L; % [rad]

% Perform least squares for each
[least_squares_ext, S_ext_2] = polyfit(T(1:index,1),gamma(1:index,1),1);
[least_squares_otw, S_otw] = polyfit(T(1:index,1),shear_strain_otw(1:index,1),1);
GJ_ext_estimate = -t/least_squares_ext(1);
GJ_otw_estimate = -t/least_squares_otw(1);

% Obtain the estimated error for each estimate
[~,ext_error_2] = polyval(least_squares_ext,T(1:index,1),S_ext_2);
[~,otw_error] = polyval(least_squares_otw,T(1:index,1),S_otw);

% Calculate theoretical GJ's
fprintf('b/t is %f therefore we can assume beta is alpha is 1/3\n\n',2*pi*R_mean/t)
otw_GJ = G*((1/3)*2*pi*R_mean*t^3);

% Calculate percent the various percent differences
perc_diff_3 = (abs(otw_GJ - GJ_ext_estimate)/GJ_ext_estimate) * 100;
fprintf('The percent difference between the estimated rigidity using the extensometer data and the rigidy using otw theory is %f\n',perc_diff_3);
perc_diff_4 = (abs(otw_GJ - GJ_otw_estimate)/GJ_otw_estimate) * 100;
fprintf('The percent difference between the estimated rigidity using the twist angle data and the rigidy using otw theory is %f\n\n',perc_diff_4);

% Create the necessary plots
figure(2)
plot(smooth(T),smooth(gamma))
hold on
plot(T,shear_strain_otw)
legend("$Extensometer\:Data$","$Twist\:Angle\:Based$",'Interpreter','latex')
xlabel("$Torque\:[lbf-in]$",'Interpreter','latex','FontSize',26)
ylabel("$Shear\:Strain[\frac{in}{in}]$",'Interpreter','latex','FontSize',26)
title("$OTW\:Section\:Shear\:Strain\:vs\:Torque$",'Interpreter','latex','FontSize',26)