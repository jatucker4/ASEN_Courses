%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
%Global formatting commands to imporve graphing looks:
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex'); 
%% Begin data processing
% Read in the text file
data = readmatrix("5mincenter_renamed");
data = data(2:end,:);
%% Next I'll get the freq spectrum for each channel's acceleration
time_vec = data(:,1) - data(1,1);
fuselage_accel = data(:,2);
tail_accel = data(:,3);
nose_accel = data(:,4);
wing_accel = data(:,5);
vibrometer_disp = data(:,10);

fuselage_disp = data(:,6);
tail_disp = data(:,7);
nose_disp = data(:,8);
wing_disp = data(:,9);

freq = interp1([data(1,1) data(end,1)], [2 50], data(1:end,1));
%% Plot the displacement and acceleration vs freq plots
figure
plot(freq,fuselage_accel)
title('Shaker Acceleration vs Frequency')
xlabel('Frequency [Hz]')
ylabel('Acceleration [$\frac{mm}{s^2}$]')

figure
plot(freq,tail_accel)
title('Tail Acceleration vs Frequency')
xlabel('Frequency [Hz]')
ylabel('Acceleration [$\frac{mm}{s^2}$]')

figure
plot(freq,nose_accel)
title('Nose Acceleration vs Frequency')
xlabel('Frequency [Hz]')
ylabel('Acceleration [$\frac{mm}{s^2}$]')

figure
plot(freq,wing_accel)
title('Wing Acceleration vs Frequency')
xlabel('Frequency [Hz]')
ylabel('Acceleration [$\frac{mm}{s^2}$]')

figure
plot(freq,vibrometer_disp)
title('Laser Vibrometer Displacement vs Frequency')
xlabel('Frequency [Hz]')
ylabel('Displacement [mm]')
%% Now perform mag factor calculations
figure
plot(freq,tail_accel./(fuselage_accel))
title('Tail Magnification Factor vs Frequency')
xlabel('Frequency [Hz]')
ylabel('Acceleration [$\frac{mm}{s^2}$]')

figure
plot(freq,nose_accel./(fuselage_accel))
title('Nose Magnification Factor vs Frequency')
xlabel('Frequency [Hz]')
ylabel('Acceleration [$\frac{mm}{s^2}$]')

figure
plot(freq,wing_accel./(fuselage_accel))
title('Wing Magnification Factor vs Frequency')
xlabel('Frequency [Hz]')
ylabel('Acceleration [$\frac{mm}{s^2}$]')

%% FFT Plots
n = 2^nextpow2(length(fuselage_accel));
Fs = freq(2);

fuselage_fft = abs(fft(fuselage_accel)/length(fuselage_accel))';
fuselage_fft = fuselage_fft(:,1:n/2+1);
fuselage_fft(:,2:end-1) = 2*fuselage_fft(:,2:end-1);

nose_fft = abs(fft(nose_accel)/length(fuselage_accel))';
nose_fft = nose_fft(:,1:n/2+1);
nose_fft(:,2:end-1) = 2*nose_fft(:,2:end-1);

tail_fft = abs(fft(tail_accel)/length(fuselage_accel))';
tail_fft = tail_fft(:,1:n/2+1);
tail_fft(:,2:end-1) = 2*tail_fft(:,2:end-1);

wing_fft = abs(fft(wing_accel)/length(fuselage_accel))';
wing_fft = wing_fft(:,1:n/2+1);
wing_fft(:,2:end-1) = 2*wing_fft(:,2:end-1);

plot_vec = linspace(2,50,length(wing_fft));
figure
plot(plot_vec,fuselage_fft)
xlabel('Frequency [Hz]')
ylabel('Amplitude')
title('System FFT Result')

% figure
hold on
plot(plot_vec,nose_fft)
% xlabel('Frequency [Hz]')
% ylabel('Amplitude')
% title('Nose FFT Result')

% figure
hold on
plot(plot_vec,tail_fft)
% xlabel('Frequency [Hz]')
% ylabel('Amplitude')
% title('Tail FFT Result')

% figure
hold on
plot(plot_vec,wing_fft)
% xlabel('Frequency [Hz]')
% ylabel('Amplitude')
% title('Wing FFT Result')