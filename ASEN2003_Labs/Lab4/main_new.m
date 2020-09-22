%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003 Lab 3: Balanced and Unbalanced Wheel
%
% Purpose: 
%           The purpose of this script is to act as the driver for the
%           model functions that calculate the angular velocity of the
%           cylinder apparatus. In addition to this it creates all
%           necessary plots.
% Authors:
%   Emma Tomlinson
%   Marlin Jacobson
%   Johnathan Tucker
%   
%
% Date Created: 02/27/2019
%
% Last Modified: 03/12/2019
%
% Inputs:
%           balanced_1: This is the data file for the first balanced trial
%           balanced_2: This is the data file for the second balanced trial
%           unbalanced_1: This is the data file for the first unbalanced
%                         trial
%           unbalanced_2: This is the data file for the second unbalanced
%                         trial
%
% Outputs:
%           This script outputs a variety of plots including the models
%           overlayed on the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clc;
clear;
close all;
%% Data Processing
% Get the data from the files
[time1,theta1,omega1] = dataProcess('balanced_1');
[time2,theta2,omega2] = dataProcess('balanced_2');
[time3,theta3,omega3] = dataProcess('unbalanced_1');
[time4,theta4,omega4] = dataProcess('unbalanced_2');
% Remove noise from data set
theta1(18) = [];
omega1(17) = 3.85;
omega1(18) = [];
time1(18) = [];
%%  Define Constants
M = 11.7; % mass of cylinder (kg)
M0 = 0.7; % mass of trailing supports (kg)
m = 3.4; % mass of extra mass (kg)
R = 0.235; % radius of cylinder (m)
k = 0.203; % radius of gyration of wheel (m)
I = M*k^2; % moment of Inertia
beta = 5.14*(pi/180); % slope of ramp (rad)
r = 0.178; % radius to extra mass (m)
rem = 0.019; % radius of extra mass (m)
g = 9.81; % gravitational acceleration (m/s^2)
h1 = 0.37338 + R; % start height of of cyclinder (m)
%% Model Verification/Plot of expected values
% Create the test theta vectory
theta = linspace(0,15);
% Calculate the model test omega vectors
omega_test_1 = model_1(theta);
omega_test_2 = model_2(0.6,theta);
omega_test_3 = MODEL_3(0.6,theta);
omega_test_4 = model_4(0.6,theta);

% Plot the test results
figure(1)
plot(theta,omega_test_1,theta,omega_test_2,theta,omega_test_3,theta,...
    omega_test_4,'LineWidth', 2)
title('$Modeled\:Angular\:Velocity\:Versus\:Angular\:Position$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(radians)$','Interpreter','latex',...
    'FontSize',26)
ylabel('$Angular\:Velocity(rad/s)$','Interpreter','latex','FontSize',26)
legend('Model 1','Model 2','Model 3', 'Model 4')
%% Begin omega calculations
% Get omega using the model 1 derivation
omega_model_1 = model_1(theta1);
% Get the angular acceleration for the model 1 calculation
alpha_model_1 = omega_model_1./theta1;
moment_1 = getMoment(omega1,theta1);

% Get omega using the model 2 derivation
omega_model_2 = model_2(moment_1,theta1);
% Get the angular acceleration for the model 2 calculation
alpha_model_2 = omega_model_2./theta1;

% Get omega using the model 3 derivation
omega_model_3 = MODEL_3(moment_1,theta3);
% Get the angular acceleration for the model 3 calculation
alpha_model_3 = omega_model_3./theta3;

% Get omega using the model 4 derivation
omega_model_4 = model_4(moment_1,theta3);
% Get the angular acceleration for the model 4 calculation
alpha_model_4 = omega_model_4./theta3;

%% Plot results

% plot balanced wheel data
figure(2);
plot(theta1,omega1,'o');
hold on;
plot(theta2,omega2,'x');
hold on
plot(theta1,omega_model_1,'--','LineWidth', 2);
hold on
plot(theta1,omega_model_2,'LineWidth', 2);

xlabel('$Angular\:Position\:Theta\:(rad)$','Interpreter','Latex',...
    'FontSize',26);
ylabel('$Angular\:Velocity\:Omega\:(rad/s)$','Interpreter','Latex',...
    'FontSize',26);
title('$Angular\:Velocity\:vs\:Angular Position\:(Balanced Wheel)$',...
    'Interpreter','Latex','FontSize',26);
legend('Trial 1','Trial 2','Model 1','Model 2');
hold off;

% plot unbalanced data

figure(3);
plot(theta3,omega3,'x');
hold on;
plot(theta4,omega4,'o');
hold on
plot(theta3,omega_model_3,'--','LineWidth', 2);
hold on
plot(theta3,omega_model_4,'LineWidth', 2)
xlabel('$Angular\:Position\:Theta\:(rad)$','Interpreter','Latex',...
    'FontSize',26);
ylabel('$Angular\:Velocity\:Omega\:(rad/s)$','Interpreter','Latex',...
    'FontSize',26);
title('$Angular\:Velocity\:vs\:Angular\:Position\:(Unbalanced Wheel)$',...
    'Interpreter','Latex','FontSize',26);
legend('Trial 1','Trial 2','Model 3','Model 4');
hold off;
%% Plot the risiduals
% First calculate them
% Get the risiduals for every model and every balanced trial
risiduals_mod1_trial1 = omega1 - omega_model_1;
risiduals_mod1_trial2 = omega2 - model_1(theta2);
risiduals_mod2_trial1 = omega1 - omega_model_2;
risiduals_mod2_trial2 = omega2 - model_2(moment_1,theta2);

% Get the risiduals for every model and every unbalanced trial
risiduals_mod3_trial1 = omega3 - omega_model_3;
risiduals_mod3_trial2 = omega4 - MODEL_3(moment_1,theta4);
risiduals_mod4_trial1 = omega3 - omega_model_4;
risiduals_mod4_trial2 = omega4 - model_4(moment_1,theta4);

figure(4)
plot(theta1, risiduals_mod1_trial1,'LineWidth', 2)
hold on
plot(theta2, risiduals_mod1_trial2,'LineWidth', 2)
xlabel('$Angular\:Position:Theta\:(rad)$','Interpreter','Latex',...
    'FontSize',26);
ylabel('$Angular\:Velocity\:Omega\:(rad/s)$','Interpreter',...
    'Latex','FontSize',26);
title('$Model\:One\:Risiduals\:(Balanced\:Wheel)$',...
    'Interpreter','Latex','FontSize',26);
legend('Trial 1 Difference','Trial 2 Difference');
hold off

figure(5)
plot(theta1, risiduals_mod2_trial1,'LineWidth', 2)
hold on
plot(theta2, risiduals_mod2_trial2,'LineWidth', 2)
xlabel('$Angular\:Position\:Theta\:(rad)$','Interpreter','Latex',...
    'FontSize',26);
ylabel('$Angular\:Velocity\:Omega\:(rad/s)$','Interpreter','Latex',...
'FontSize',26);
title('$Model\:Two\:Risiduals\:(Balanced\:Wheel)$',...
    'Interpreter','Latex','FontSize',26);
legend('Trial 1 Difference','Trial 2 Difference');
hold off

figure(6)
plot(theta3, risiduals_mod3_trial1,'LineWidth', 2)
hold on
plot(theta4, risiduals_mod3_trial2,'LineWidth', 2)
xlabel('$Angular\:Position\:Theta\:(rad)$','Interpreter','Latex',...
    'FontSize',26);
ylabel('$Angular\:Velocity\:Omega\:(rad/s)$','Interpreter','Latex',...
'FontSize',26);
title('$Model\:Three\:Risiduals\:(Unbalanced\:Wheel)$',...
    'Interpreter','Latex','FontSize',26);
legend('Trial 1 Difference','Trial 2 Difference');
hold off

figure(7)
plot(theta3, risiduals_mod4_trial1,'LineWidth', 2)
hold on
plot(theta4, risiduals_mod4_trial2,'LineWidth', 2)
xlabel('$Angular\:Position\:Theta\:(rad)$','Interpreter','Latex',...
    'FontSize',26);
ylabel('$Angular\:Velocity\:Omega\:(rad/s)$','Interpreter','Latex',...
'FontSize',26);
title('$Model\:Four\:Risiduals\:(Unalanced\:Wheel)$',...
    'Interpreter','Latex','FontSize',26);
legend('Trial 1 Difference','Trial 2 Difference');
hold off
%% Calculate Statistics
% Get the statistics for every model and balanaced trial
statistics_mod1_trial1 = getStatistics(risiduals_mod1_trial1);
statistics_mod1_trial2 = getStatistics(risiduals_mod1_trial2);
statistics_mod2_trial1 = getStatistics(risiduals_mod2_trial1);
statistics_mod2_trial2 = getStatistics(risiduals_mod2_trial2);

% Get the statisitcs for every model and unbalanced trial
statistics_mod3_trial1 = getStatistics(risiduals_mod3_trial1);
statistics_mod3_trial2 = getStatistics(risiduals_mod3_trial2);
statistics_mod4_trial1 = getStatistics(risiduals_mod4_trial1);
statistics_mod4_trial2 = getStatistics(risiduals_mod4_trial2);

%% Plot the angular acceleration
figure(8)
plot(theta1,alpha_model_1,'LineWidth', 2)
hold on
plot(theta1,alpha_model_2,'LineWidth', 2)
hold on
plot(theta3,alpha_model_3,'LineWidth', 2)
hold on
plot(theta3,alpha_model_3,'LineWidth', 2)
xlabel('$Angular\:Position\:Theta\:(rad)$','Interpreter','Latex',...
    'FontSize',26);
ylabel('$Angular\:Acceleration\:Alpha\:(rad/s^2)$','Interpreter','Latex',...
'FontSize',26);
title('$Angular\:Acceleration\:vs\:Angular\:Position$',...
    'Interpreter','Latex','FontSize',26);
legend('Model 1','Model 2','Model 3','Model 4');
 
alpha_model_1_1 = mean(alpha_model_1);
alpha_model_1_2 = mean(model_1(theta2)./theta2);
alpha_model_2_1 = mean(alpha_model_2);
alpha_model_2_2 = mean(model_2(moment_1,theta2)./theta2);

% Get the constant accelerations
accel_model_1_1 = alpha_model_1_1*R;
accel_model_1_2 = alpha_model_1_2*R;
accel_model_2_1 = alpha_model_2_1*R;
accel_model_2_2 = alpha_model_2_2*R;
%% Get the table of moment values
moment_1 = model_2(1.0,theta1);
moment_2 = model_2(0.9,theta1);
moment_3 = model_2(0.8,theta1);
moment_4 = model_2(0.7,theta1);
moment_5 = model_2(0.6,theta1);

figure(9)
plot(theta1,moment_1,theta1,moment_2,theta1,moment_3,theta1,moment_4,...
    theta1,moment_5,theta1,omega_model_2,'LineWidth', 2)
xlabel('$Angular\:Position\:Theta\:(rad)$','Interpreter','Latex',...
    'FontSize',26);
ylabel('$Angular\:Velocity\:Omega\:(rad/s)$','Interpreter','Latex',...
'FontSize',26);
title('$Angular\:Velocity\:Under\:Varying\:Moments$',...
    'Interpreter','Latex','FontSize',26);
legend('Moment = 1.0 N-m','Moment = 0.9 N-m','Moment = 0.8 N-m',...
    'Moment = 0.7 N-m','Moment = 0.6 N-m','Experimental Data');






