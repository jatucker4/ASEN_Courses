%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003 Lab 3: Balanced and Unbalanced Wheel

% Purpose: 

% Authors:
%   Emma Tomlinson
%   
%   
%   

% Date Created: 02/27/2019

% Last Modified: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% housekeeping
clc;
clear all;
close all;

% Model
theta = [0:15];
modOmega = modelAngVel(theta);

% Experiment
filename1 = 'balanced_1.dms';
filename2 = 'balanced_2.dms';
filename3 = 'unbalanced_1.dms';
filename4 = 'unbalanced_2.dms';
%alpha = expAngVel(filename1,filename2,filename3,filename4);
[alpha1,theta1,omega1] = angAccel(filename1);
[alpha2,theta2,omega2] = angAccel(filename2);
[alpha3,theta3,omega3] = angAccel(filename3);
[alpha4,theta4,omega4] = angAccel(filename4);

% plot balanced wheel data
figure;
plot(theta1,omega1,'--');
hold on;
plot(theta2,omega2);
xlabel('Angular Position Theta (rad)');
ylabel('Angular Velocity Omega (rad/s)');
title('Angular Velocity vs Angular Position (Balanced Wheel)');
legend('Trial 1','Trial 2');
hold off;

% plot unbalanced data

figure;
plot(theta3,omega3,'--');
hold on;
plot(theta4,omega4);
xlabel('Angular Position Theta (rad)');
ylabel('Angular Velocity Omega (rad/s)');
title('Angular Velocity vs Angular Position (Unbalanced Wheel)');
legend('Trial 1','Trial 2');
hold off;

















