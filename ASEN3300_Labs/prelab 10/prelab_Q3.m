%% Housekeeping
clc;
close all; 
clear all;
%Global formatting commands to imporve graphing looks:
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex'); 
%% Create Constants
R = 33e3; 
% R = 3.3e3;
C1 = 0.01e-6;
% C1 = 1e-12;
C2 = 100e-12;
Q = 0.5*sqrt(C1/C2);
omega_0 = 1/(R*sqrt(C1*C2));
freq = 0:1:60e3;
omega = 2*pi*freq;
%% Get the Transfer Function Magnitude and Phase
transfer = (omega_0^2)./((omega_0^2 -omega.^2) + (1i*omega_0*omega/Q));
mag = abs(transfer);
phase = (180/pi)*angle(transfer);

%% Create all plots
figure
plot(freq,20*log10(mag))
hold on
ylabel('Magnitude [dBV]')
xlabel('Frequency [Hz]')
title('Filter Magnitude Response');
grid
axis tight


%% Create Constants
R = 82.5e3; 
% R = 3.3e3;
C1 = 0.016e-6;
% C1 = 1e-12;
C2 = 10e-12;
Q = 0.5*sqrt(C1/C2);
omega_0 = 1/(R*sqrt(C1*C2));
freq = 0:1:60e3;
omega = 2*pi*freq;
%% Get the Transfer Function Magnitude and Phase
transfer = (omega_0^2)./((omega_0^2 -omega.^2) + (1i*omega_0*omega/Q));
mag = abs(transfer);
phase = (180/pi)*angle(transfer);

%% Create all plots

plot(freq,20*log10(mag))
% ylabel('Magnitude [dB]')
% xlabel('Frequency [Hz]')
% title('Filter Magnitude Response');
% grid
% axis tight
legend('Q=5','Q=20')
figure
plot(freq,phase)
ylabel('Phase [deg]')
xlabel('Frequency [Hz]')
title('Filter Phase Response')
grid
axis tight

freq_vec = [1000,4000,4898,6000,7000,8000];
mag_vec = [1,3.38,4.54,2,1.2,.7];
mag_vec = 20*log10(mag_vec);
figure
plot(freq,20*log10(mag))
hold on
scatter(freq_vec,mag_vec)
ylabel('Magnitude [dBV]')
xlabel('Frequency [Hz]')
title('Filter Magnitude Response');
legend('Simulated Response','Experiment Data')
grid
axis tight

