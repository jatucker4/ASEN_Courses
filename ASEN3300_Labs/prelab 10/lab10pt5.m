clear; clc; close all;

load('asen3300mod.mat');
%signal = signalnoisy;
period = 1/fs;

x = linspace(0,period*9.75,40);

figure
plot(x,signal(1:40));
xlim([0 x(40)])
title('Input Signal vs Time')
ylabel('Input Signal')
xlabel('Time (seconds)')

FFT_sig = fft(signal);
L = length(signal);
P2 = abs(FFT_sig/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
figure
plot(f,P1)
title('Frequency Spectrum of Input Signal')
xlabel('Frequency (Hz)')
ylabel('|P1(f)|')

% figure();
% plot(x,signal(1:40));
% xlim([0 x(40)])
% title('Signal vs Time')
% ylabel('Signal')
% xlabel('Time (seconds)')


am = demod(signal,fc,fs,'am');
fm = demod(signal,fc,fs,'fm');

figure
plot(x,am(1:40));
xlim([0 x(40)])
title('Signal vs Time AM Demodulation')
ylabel('Signal')
xlabel('Time (seconds)')

figure
plot(x,fm(1:40));
xlim([0 x(40)])
title('Signal vs Time FM Demodulation')
ylabel('Signal')
xlabel('Time (seconds)')


% sound(am,fs) %This sounds bad
sound(fm,fs);

figure

Y = fft(fm);
L = length(fm);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = fs*(0:(L/2))/L;
plot(f,P1) 
xlim([0 5000])
title('Frequency Spectrum After FM Demodulation')
xlabel('Frequency (Hz)')
ylabel('|P1(f)|')


