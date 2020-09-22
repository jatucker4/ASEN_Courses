%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script was created for the fourth lab in ASEN 3300. It loads in
% three different signal mat files and performs frequency analysis on them
% using the fft function.
% 
% Created by: Johnathan Tucker
%
% Created on: 2/9/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
close all;
clear all;
%% Perform analysis on the first signal
load('signal1.mat')

% Plot the first signal
figure(1)
subplot(1,2,1)
plot(t(1:100),x(1:100))
xlabel('$Time\:[s]$','Interpreter','latex','FontSize',26)
ylabel('$Amplitude\:[V]$','Interpreter','latex','FontSize',26)
sgtitle('$Signal\:And\:Corresponding\:Frequency\:Spectrum$',...
    'Interpreter','latex','FontSize',26)
% Compute Fourier Transform of the signal 
FFT_1 = fft(x);
FFT_1 = abs(FFT_1/length(FFT_1));
f = Fs*(0:length(x)-1)/length(x);
subplot(1,2,2)
plot(f,FFT_1)
xlabel('$Frequency\:[Hz]$','Interpreter','latex','FontSize',26)
ylabel('$Amplitude\:[dBV_{RMS}]$','Interpreter','latex','FontSize',26)
%% Perform analysis on the second signal
load('signal2.mat')

% Plot the first signal
figure(2)
subplot(1,2,1)
plot(t(1:100),x(1:100))
xlabel('$Time\:[s]$','Interpreter','latex','FontSize',26)
ylabel('$Amplitude\:[V]$','Interpreter','latex','FontSize',26)
sgtitle('$Signal\:And\:Corresponding\:Frequency\:Spectrum$',...
    'Interpreter','latex','FontSize',26)

% Compute Fourier Transform of the signal 
FFT_1 = fft(x);
FFT_1 = abs(FFT_1/length(FFT_1));
f = Fs*(0:length(x)-1)/length(x);
subplot(1,2,2)
plot(f,FFT_1)
xlabel('$Frequency\:[Hz]$','Interpreter','latex','FontSize',26)
ylabel('$Amplitude\:[dBV_{RMS}]$','Interpreter','latex','FontSize',26)
%% Perform analysis on the third signal
load('signal3.mat')

% Plot the first signal
figure(3)
subplot(1,2,1)
plot(t(1:100),x(1:100))
xlabel('$Time\:[s]$','Interpreter','latex','FontSize',26)
ylabel('$Amplitude\:[V]$','Interpreter','latex','FontSize',26)
sgtitle('$Signal\:And\:Corresponding\:Frequency\:Spectrum$',...
    'Interpreter','latex','FontSize',26)

% Compute Fourier Transform of the signal 
FFT_1 = fft(x);
FFT_1 = abs(FFT_1/length(FFT_1));
f = Fs*(0:length(x)-1)/length(x);
subplot(1,2,2)
plot(f,FFT_1)
xlabel('$Frequency\:[Hz]$','Interpreter','latex','FontSize',26)
ylabel('$Amplitude\:[dBV_{RMS}]$','Interpreter','latex','FontSize',26)