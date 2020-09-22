% ASEN 2003 - Lab 5 - YoYo Despinner
% 3/18/2019 - Group 2 - Jashan Chopra, John Tucker, Davis Peirce, Andrew Fu
% This script plots the derivation and actual data for the second lab over
% time in order to compare the predicited and actual data

clc; close all; clear; 

%% Define Initial Variables
w0 = 130; % initial angular velocity [rpm]
w0 = (w0 * 2 * pi) / 60; % Convert to [rad/s]

I = .0063; % Moment of Inertia -- Note Uncertainity for later
rOuter = .0762; % Outer Radius [m]
m1 = 54 / 1000; m2 = m1; % Masses [kg] 

C = (I / ((m1+m2)*rOuter^2)) + 1; % Variable C (Moment of Inertia Normalization?)
totTime = sqrt(C/(w0^2)); % Total time to stop
%% Length of Cord & Total Time [Tangential Despin Prediction]
totLength = sqrt(C*rOuter^2); % Total length of the cord
length = linspace(0,totLength,1000); % created length vector (totLength is the calc final length)

top2 = (C*rOuter^2) - length.^2;
bot2 = (C*rOuter^2) + length.^2;
wLength = w0*(top2./bot2); % [rad/s] - angular velocity (func of length)

%% Length of Cord [Radial Despin]

totLengthRadial = rOuter*(sqrt(C)-1); % Total length of cord radially found
lengthRadial = linspace(0,totLengthRadial,1000); % created lengthRadial

% ! THIS FUNCTION NEEDS TO BE DIFFERENT ! % 
topRad = (C*rOuter^2) - lengthRadial.^2;
botRad = (C*rOuter^2) + lengthRadial.^2;
wLengthRadial = (-w0*(I+m1*rOuter^2)+(rOuter*w0*sqrt(C)*...
    (lengthRadial+rOuter)))/I; % [rad/s] - angular velocity (func of length)

%% Derived Equations[Tangential Despin]

time = linspace(0,totTime,1000); % created time vector (totTime is the calculated final time) 

% Derivation 1 - Angular Velocity
    
top1 = C - (w0^2*time.^2); % Numerator of derivation
bot1 = C + (w0^2*time.^2); % Denominator of derivation
wTime = w0*(top1./bot1); % [rad/s] - angular velocity (func of time)

% Derivation 2 - Angular Acceleration
top3 = (-4 * w0^2 * time * C);
bot3 = (C + w0^2 * time.^2).^2;
alpha = w0*(top3./bot3); % [rad/s^2] - angular acceleration (func of time) 

% Derivation 3 - Cable Tension
tension = abs((I / (2*rOuter)) * alpha); % absolute value of cable tension (func of alpha)

%% Plots [Tangential Derivation]

figure(1)
plot(time,wTime,'LineWidth',2)
% Angular velocity plot
title('Angular Velocity Versus Time',...
    'FontSize',26); 
xlabel('Time[s]','FontSize',26);
ylabel('Angular Velocity[rad/s]','FontSize',26);
legend('Angular Velocity-Tangential Release');

figure(2)
plot(time,alpha) % Angular Accel Plot
title('Angular Acceleration Versus Time','FontSize',26); 
ylabel('Angular Acceleration [rad/s^2]','FontSize',26)
xlabel('Time[s]','FontSize',26);
legend('Angular Acceleration-Tangential Release');

figure(3)
plot(time,tension) % Cable Tension Plot
title('Cable Tension over Time','FontSize',26); 
xlabel('Time [s]','FontSize',26); 
ylabel('Cable Tension [N]','FontSize',26)
legend('Tension-Tangential Release');
%% Calculate the Moment due to friction

%% Plot [Angular Velocity v length]
figure(4)
plot(length,wLength) % Tangential (func of length)

hold on
plot(lengthRadial,wLengthRadial)% Radial (func of length)

title('Angular Velocity v. Length')
xlabel('Length [m]')
ylabel('Angular Velocity [rad/s]')
legend('Tangential Derivation', 'Radial Derivation')

%% Importating Data 
filename1 = 0;
% data1 = xlsread(filename1);

%% Compare Residuals between Actual and Predicted Data
%% Theory section plot
