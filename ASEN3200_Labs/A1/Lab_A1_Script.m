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
%% Import automatic test data
% Test one
automatic_1 = readmatrix('automatic_rw_0_75A_0_5Hz.csv');
% Test two
automatic_2 = readmatrix('automatic_rw_0_25A_0_1Hz.csv');
% Test three
automatic_3 = readmatrix('automatic_rw_0_5A_0_2Hz.csv');

%% Convert the input to rad/s
% Test one
automatic_1(:,3) = automatic_1(:,3)*(2*pi/60);
% Test two
automatic_2(:,3) = automatic_2(:,3)*(2*pi/60);
% Test three
automatic_3(:,3) = automatic_3(:,3)*(2*pi/60);
%% Calculate the slope(scale factor) and bias (y-intercept) for each test
% Test one
p_auto_1 = polyfit(automatic_1(:,3),automatic_1(:,2),1);
% Test two
p_auto_2 = polyfit(automatic_2(:,3),automatic_2(:,2),1);
% Test three
p_auto_3 = polyfit(automatic_3(:,3),automatic_3(:,2),1);
%% Create best fit lines for each test
% Test one
test_one_line = p_auto_1(1)*automatic_1(:,3) + p_auto_1(2);

% Test two
test_two_line = p_auto_2(1)*automatic_2(:,3) + p_auto_2(2);

% Test three
test_three_line = p_auto_3(1)*automatic_3(:,3) + p_auto_3(2);
%% Create plots for each test including best fit line
% Test one
figure(1)
plot(automatic_1(:,3),automatic_1(:,2))
hold on
plot(automatic_1(:,3),test_one_line)
xlabel("$Encoder\:Velocity\:[rad/s]$",'Interpreter','latex','FontSize',26)
ylabel("$Gyro\:Output\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Rate\:Gyro\:With\:0.75A\:and\:0.5Hz\:Control$",'Interpreter','latex','FontSize',26)
legend("Experimental Data", "Best Fit Line")

% Test two
figure(2)
plot(automatic_2(:,3),automatic_2(:,2))
hold on
plot(automatic_2(:,3),test_two_line)
xlabel("$Encoder\:Velocity\:[rad/s]$",'Interpreter','latex','FontSize',26)
ylabel("$Gyro\:Output\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Rate\:Gyro\:With\:0.25A\:and\:0.1Hz\:Control$",'Interpreter','latex','FontSize',26)
legend("Experimental Data", "Best Fit Line")

% Test three
figure(3)
plot(automatic_3(:,3),automatic_3(:,2))
hold on
plot(automatic_3(:,3),test_three_line)
xlabel("$Encoder\:Velocity\:[rad/s]$",'Interpreter','latex','FontSize',26)
ylabel("$Gyro\:Output\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Rate\:Gyro\:With\:0.5A\:and\:0.2Hz\:Control$",'Interpreter','latex','FontSize',26)
legend("Experimental Data", "Best Fit Line")
%% Calculating the mean and std for bias and sensitivity
mean_sens = mean([p_auto_1(1),p_auto_2(1),p_auto_3(1)]);
std_sens = std([p_auto_1(1),p_auto_2(1),p_auto_3(1)]);
mean_bias = mean([p_auto_1(2),p_auto_2(2),p_auto_3(2)]);
std_bias = std([p_auto_1(1),p_auto_2(1),p_auto_3(1)]);
%% Ang vel over time exp versus analytical
figure(4)
subplot(2,1,1)
plot(automatic_1(1552:end,1)-automatic_1(1552,1),automatic_1(1552:end,2))
xlabel("$Time\:[s]$",'Interpreter','latex')
ylabel("$Angular\:Velocity\:[rad/s]$",'Interpreter','latex')
title("$Experimental\:Data:$",'Interpreter','latex')

subplot(2,1,2)
plot(automatic_1(1552:end,1)-automatic_1(1552,1),automatic_1(1552:end,3).*mean_sens)
xlabel("$Time\:[s]$",'Interpreter','latex')
ylabel("$Angular\:Velocity\:[rad/s]$",'Interpreter','latex')
title("$Analytical\:Data:$",'Interpreter','latex')
sgtitle("$0.74A\:0.5Hz\:Angular\:Velocity\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

figure(5)
subplot(2,1,1)
plot(automatic_2(2682:end,1)-automatic_2(2682,1),automatic_2(2682:end,2))
xlabel("$Time\:[s]$",'Interpreter','latex')
ylabel("$Angular\:Velocity\:[rad/s]$",'Interpreter','latex')
title("$Experimental\:Data:$",'Interpreter','latex')

subplot(2,1,2)
plot(automatic_2(2682:end,1)-automatic_2(2682,1),automatic_2(2682:end,3).*mean_sens)
xlabel("$Time\:[s]$",'Interpreter','latex')
ylabel("$Angular\:Velocity\:[rad/s]$",'Interpreter','latex')
title("$Analytical\:Data:$",'Interpreter','latex')
sgtitle("$0.25A\:0.1Hz\:Angular\:Velocity\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

figure(6)
subplot(2,1,1)
plot(automatic_3(2:end,1)-automatic_3(2,1),automatic_3(2:end,2))
xlabel("$Time\:[s]$",'Interpreter','latex')
ylabel("$Angular\:Velocity\:[rad/s]$",'Interpreter','latex')
title("$Experimental\:Data:$",'Interpreter','latex')

subplot(2,1,2)
plot(automatic_3(2:end,1)-automatic_3(2,1),automatic_3(2:end,3).*mean_sens)
xlabel("$Time\:[s]$",'Interpreter','latex')
ylabel("$Angular\:Velocity\:[rad/s]$",'Interpreter','latex')
title("$Analytical\:Data:$",'Interpreter','latex')
sgtitle("$0.5A\:0.2Hz\:Angular\:Velocity\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)
%% Get the angular position and plot it
delta_t_1 = automatic_1(2,1) - automatic_1(3,1);
delta_t_2 = automatic_2(2,1) - automatic_2(3,1);
delta_t_3 = automatic_3(2,1) - automatic_3(3,1);

figure(7)
plot(automatic_1(1552:end,1)-automatic_1(1552,1),automatic_1(1552:end,3).*(mean_sens*delta_t_1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',15)
ylabel("$Angular\:Position\:[rad]$",'Interpreter','latex','FontSize',15)
title("$Angular\:Position\:Over\:Time\:With\:0.75A\:and\:0.5Hz\:Control$",'Interpreter','latex','FontSize',15)

figure(8)
plot(automatic_2(2682:end,1)-automatic_2(2682,1),automatic_2(2682:end,3).*(mean_sens*delta_t_2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',15)
ylabel("$Angular\:Position\:[rad]$",'Interpreter','latex','FontSize',15)
title("$Angular\:Position\:Over\:Time\:With\:0.25A\:and\:0.1Hz\:Control$",'Interpreter','latex','FontSize',15)

figure(9)
plot(automatic_3(2:end,1)-automatic_3(2,1),automatic_3(2:end,3).*(mean_sens*delta_t_3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',15)
ylabel("$Angular\:Position\:[rad]$",'Interpreter','latex','FontSize',15)
title("$Angular\:Position\:Over\:Time\:With\:0.5A\:and\:0.2Hz\:Control$",'Interpreter','latex','FontSize',15)
%% Question 2 Import data
% Test 1
torque_1 = readmatrix('torque_1_10mNm.csv');
torque_1(:,3) = torque_1(:,3)*(2*pi/60);
% Test 2
torque_2 = readmatrix('torque_2_15mNm.csv');
torque_2(:,3) = torque_2(:,3)*(2*pi/60);
% Test 3
torque_3 = readmatrix('torque_3_20mNm.csv');
torque_3(:,3) = torque_3(:,3)*(2*pi/60);
% Test 4
torque_4 = readmatrix('torque_4_5mNm.csv');
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
I_1 = (10e-3)/(p_torque_1(1));
I_2 = (15e-3)/(p_torque_2(1));
I_3 = (20e-3)/(p_torque_3(1));
I_4 = (5e-3)/(p_torque_4(1));

mean_I = mean([I_1,I_2,I_3,I_4]);
std_I = std([I_1,I_2,I_3,I_4]);
