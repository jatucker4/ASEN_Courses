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

%{
lengthRadial = linspace(0,totLengthRadial,1000); % created lengthRadial
topRad = (C*rOuter^2) - lengthRadial.^2;
botRad = (C*rOuter^2) + lengthRadial.^2;
wLengthRadial = w0*(topRad./botRad); % [rad/s] - angular velocity (func of length)
%}

%% Derived Equations[Tangential Despin]

time = linspace(0,totTime,28); % created time vector (totTime is the calculated final time) 

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
plot(time,wTime) % Angular velocity plot [Predicted]
hold on
title('Angular Velocity over time','fontsize',24); xlabel('Time [s]','fontsize',24); 
ylabel('Angular Velocity [rad/s]','fontsize',24)

figure(2)
plot(time,alpha) % Angular Accel Plot [Predicted]
title('Angular Acceleration over time','fontsize',24); xlabel('Time [s]','fontsize',24); 
ylabel('Angular Acceleration [rad/s^2]','fontsize',24)

figure(3)
plot(time,tension) % Cable Tension Plot [Predicted]
title('Cable Tension over Time','fontsize',24); xlabel('Time [s]','fontsize',24); 
ylabel('Cable Tension [N]','fontsize',24)

%% Plot [Angular Velocity v length]
figure(4)
plot(length,wLength) % Tangential (func of length) [Predicted]
title('Angular Velocity v. Length','fontsize',24)
xlabel('Length [m]','fontsize',24)
ylabel('Angular Velocity [rad/s]','fontsize',24)

%% Importating Data 
filename1 = 'mass_trial_2'; % Filenames
filename2 = 'no_mass_trial';

dataMass = load(filename1); % Initial load
dataNoMass = load(filename2);

dataMass = dataMass(15:42,:); % Aquire correct points
dataNoMass = dataNoMass(2:3900,:);

dataMass(:,2) = dataMass(:,2) * (2*pi/60); % Convert to rads
dataNoMass(:,2) = dataNoMass(:,2) * (2*pi/60);

dataMass(:,1) = (dataMass(:,1) / 1000) - .14; % convert to seconds and fix to 0
dataNoMass(:,1) = (dataNoMass(:,1) / 1000) - .14;

figure(5) % Plot on same graph
plot(dataNoMass(:,1),dataNoMass(:,2))
hold on
plot(dataMass(:,1),dataMass(:,2))
title('Experimental Angular Velocity','fontsize',24)
legend('No Masses (Only Friction)', 'With Masses')
xlabel('Time [s]','fontsize',24)
ylabel('Angular Velocity [rad/s]','fontsize',24)

%% Compare Residuals between Actual and Predicted Data

% Calculate the angular velocity from derivation but with the same
% timescale as the experiment instead of the calculated total time
top1New = C - (w0^2*dataMass(:,1).^2);
bot1New = C + (w0^2*dataMass(:,1).^2);
wTimeNew = w0*(top1New./bot1New); % [rad/s] - angular velocity (func of time)

% Plot the experimental against the derived
figure(6)
plot(time(13),wTime(13),'o')
hold on
plot(time,wTime) 
hold on
plot(dataMass(:,1),dataMass(:,2))
title('Experimental v. Derived Angular Velocity','fontsize',24)
xlabel('Time [s]','fontsize',24)
ylabel('Angular Velocity [rad/s]','fontsize',24)
legend('Point of max length','Derived', 'Experimental')

% Plot the residuals
residual = (dataMass(:,2)' - wTime);
figure(7)
plot(dataMass(:,1),residual)
title('Angular Velocity Residuals','fontsize',24)
xlabel('Time [s]','fontsize',24)
ylabel('Angular Velocity [rad/s]','fontsize',24)

%% Compute angular acceleration and cable tension of experimental data

% Angular Acceleration [Experimental]
dwdt = gradient(dataMass(:,2)) ./ gradient(dataMass(:,1));

% Cable Tension [Experimental]
tension_ex = abs((I / (2*rOuter)) * dwdt); % absolute value of cable tension (func of alpha)

% find max points for angular acceleration
[M,I] = max(abs(alpha));
timeMaxDer = time(I);
[M,I] = max(abs(dwdt));
timeMaxExp = dataMass(I,1);

% Find max points for cable tension
[M,I] = max(abs(tension));
timeMaxDerT = time(I);
[M,I] = max(abs(tension_ex));
timeMaxExpT = dataMass(I,1);

figure(8)
plot(time,alpha)
hold on
plot(dataMass(:,1),dwdt)
plot(timeMaxDer,min((alpha)),'o')
plot(timeMaxExp,min((dwdt)),'o')
title('Angular Acceleration Comparison','fontsize',24)
xlabel('Time [s]','fontsize',24)
ylabel('Angular Acceleration [rad/s^2]','fontsize',24)
legend('Derived','Experimental','Max Derived', 'Max Eperimental')

figure(9)
plot(time,tension)
hold on
plot(dataMass(:,1),tension_ex)
plot(timeMaxDerT,max(abs(tension)),'o')
plot(timeMaxExpT,max(abs(tension_ex)),'o')
title('Cable Tension Comparison','fontsize',24)
xlabel('Time [s]','fontsize',24)
ylabel('Cable Tension [N]','fontsize',24)
legend('Derived','Experimental','Max Derived', 'Max Eperimental')