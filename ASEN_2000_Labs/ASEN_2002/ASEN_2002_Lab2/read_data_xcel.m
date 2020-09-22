function [atm_pressure,atm_temp, airspeed, auxillary, eld_x, eld_y, Voltage] = read_data_xcel(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%[atm_pressure,atm_temp, airspeed, auxillary, eld_x, eld_y, Voltage] = xlsread(filename);
data = csvread(filename, 1, 0);
atm_pressure = data(:,1);
atm_temp = data(:,2);
airspeed = data(:,3);
auxillary = data(:,4);
eld_x = data(:,5);
eld_y = data(:,6);
Voltage = data(:,7);
end

