%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver script for the first attitude lab in ASEN 3200.
%
% Date created: 11/13/2019
%
% Created by: John Tucker
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;

%% Question 2 Import data
% Test 1
torque_1 = readmatrix('Task1T12D2.csv');
torque_1(:,3) = torque_1(:,3)*(2*pi/60);
% Test 2
torque_2 = readmatrix('Task1T18D8.csv');
torque_2(:,3) = torque_2(:,3)*(2*pi/60);
% Test 3
torque_3 = readmatrix('Task1T16D6.csv');
torque_3(:,3) = torque_3(:,3)*(2*pi/60);
% Test 4
torque_4 = readmatrix('Task1T15D5.csv');
torque_4(:,3) = torque_4(:,3)*(2*pi/60);

%% Create a line of best fit
p_torque_1 = polyfit(torque_1(115:621,1),torque_1(115:621,3),1);
p_line_1 = p_torque_1(1)*torque_1(115:621,1) + p_torque_1(2);

p_torque_2 = polyfit(torque_2(120:528,1),torque_2(120:528,3),1);
p_line_2 = p_torque_2(1)*torque_2(120:528,1) + p_torque_2(2);

p_torque_3 = polyfit(torque_3(113:386,1),torque_3(113:386,3),1);
p_line_3 = p_torque_3(1)*torque_3(113:386,1) + p_torque_3(2);

p_torque_4 = polyfit(torque_4(119:618,1),torque_4(119:618,3),1);
p_line_4 = p_torque_4(1)*torque_4(119:618,1) + p_torque_4(2);

%% Create plots 
figure(10)
plot(torque_1(115:621,1)-torque_1(115,1),torque_1(115:621,3))
hold on
plot(torque_1(115:621,1)-torque_1(115,1),p_line_1)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Angular\:Velocity\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Angular\:Velocity\:Over\:Time\:with\:10mNm\:Torque$",'Interpreter','latex','FontSize',26)
legend("Experimental Data", "Best Fit Line")

figure(11)
plot(torque_2(120:528,1)-torque_2(120,1),torque_2(120:528,3))
hold on
plot(torque_2(120:528,1)-torque_2(120,1),p_line_2)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Angular\:Velocity\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Angular\:Velocity\:Over\:Time\:with\:15mNm\:Torque$",'Interpreter','latex','FontSize',26)
legend("Experimental Data", "Best Fit Line")

figure(12)
plot(torque_3(113:386,1)-torque_3(113,1),torque_3(113:386,3))
hold on
plot(torque_3(113:386,1)-torque_3(113,1),p_line_3)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Angular\:Velocity\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Angular\:Velocity\:Over\:Time\:with\:20mNm\:Torque$",'Interpreter','latex','FontSize',26)
legend("Experimental Data", "Best Fit Line")

figure(13)
plot(torque_4(119:618,1)-torque_4(119,1),torque_4(119:618,3))
hold on
plot(torque_4(119:618,1)-torque_4(119,1),p_line_4)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Angular\:Velocity\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Angular\:Velocity\:Over\:Time\:with\:5mNm\:Torque$",'Interpreter','latex','FontSize',26)
legend("Experimental Data", "Best Fit Line")

%% Compute moment of inertia as well as mean and std
I_1 = (12e-3)/(p_torque_1(1));
I_2 = (8e-3)/(p_torque_2(1));
I_3 = (6e-3)/(p_torque_3(1));
I_4 = (5e-3)/(p_torque_4(1));

mean_I = mean([I_1,I_2,I_3,I_4]);
std_I = std([I_1,I_2,I_3,I_4]);

%% Using the mean I calculate the gains
zeta = 0.7329;
omega_n = 2.728;

K1 = omega_n^2 * mean_I;
K2 = omega_n * 2 * zeta * mean_I;
