clear all
clc
close all

MAX_VOLTAGE = 3.3;
MIN_VOLTAGE = 0;

voltages = [0:0.25:3.25]; % voltage levels to convert to digital signals
bits     = [4, 8, 12];    % ADC bits

for b = bits
  figure; hold on; grid on;
  title(sprintf('ADC bin vs voltage - %d bits', b));
  xlabel('voltage')
  ylabel('bin')

  for v = voltages
    bin = Voltage2Bin(v, MAX_VOLTAGE, MIN_VOLTAGE, b);
    plot(v, bin, 'ro');
  end
end

% sine wave
offset      = 1.65;
amplitude   = 3.3/2;
sample_freq = 20;
delta_t     = 1/sample_freq;
t0          = 0;
tmax        = 1; % 1 period
time        = [t0:delta_t:tmax];
voltage     = offset + amplitude*sin(2*pi*time);

figure; hold on; grid on;
title(sprintf('ADC reading vs sine wave, f = %d Hz, f_s = %d Hz', 1, sample_freq));
xlabel('reaAding number');
ylabel('bin number');

for i = 1:length(voltage)
  v = voltage(i);
  % 12 bit, same logic level as before?
  bin = Voltage2Bin(v, MAX_VOLTAGE, MIN_VOLTAGE, 12); 
  plot(i, bin, 'ro');
end