%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3300 - Lab 12
% 
% Created By: Johnathan Tucker
%
% Collaborators: 
%
% The purpose of the script is to act as a driver that will execute the
% functions necesary to solve questions for ASEN 3300 Lab 12.
%
% Created Date: 4/28/2020
%
% Change Log: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
tic
%Global formatting commands to imporve graphing looks:
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex'); 
%% Begin code block to create plots
lat_park_1 = [39.830563, 39.830585, 39.830587, 39.830578, 39.830578]; 
lon_park_1 = [-105.048455, -105.048389, -105.048393, -105.048372, -105.048368];
lat_park_2 = [39.830563,39.830567, 39.830554, 39.830556, 39.830572];
lon_park_2 = [-105.048413, -105.048488, -105.048431, -105.048354, -105.048356];

lat_house_1 = [39.830452,39.830447,39.830431,39.830375,39.830369];
lon_house_1 = [-105.046155,-105.046127,-105.046145,-105.046093,-105.046110];
lat_house_2 = [39.830369,39.830381,39.830453,39.830402,39.830422];
lon_house_2 = [-105.046043, -105.045919, -105.046124, -105.046098, -105.046116];

% Average the lat long for the data points
lat_park_avg = mean([lat_park_1,lat_park_2]);
lat_park_std = std([lat_park_1,lat_park_2]);
lon_park_avg = mean([lon_park_1,lon_park_2]);
lon_park_std = std([lon_park_1,lon_park_2]);

lat_house_avg = mean([lat_house_1,lat_house_2]);
lat_house_std = std([lat_house_1,lat_house_2]);
lon_house_avg = mean([lon_house_1,lon_house_2]);
lon_house_std = std([lon_house_1,lon_house_2]);
% True locations from google maps
lat_park_true = 39.83056; 
lon_park_true = -105.048371;
lat_house_true = 39.830431;
lon_house_true = -105.046122;
dist_park = lldistkm([lat_park_avg,lon_park_avg],[lat_park_true,lon_park_true]);
dist_house = lldistkm([lat_house_avg,lon_house_avg],[lat_house_true,lon_house_true]);

figure
hold on
plot(lon_park_1,lat_park_1,'.r','MarkerSize',20)
plot(lon_park_2,lat_park_2,'.b','MarkerSize',20)
plot(lon_park_true,lat_park_true,'*k','MarkerSize',20)
hold off
xlabel('Longitude [$^\circ$]')
ylabel('Latitude [$^\circ$]')
title('GPS Application Measurements and Google Maps Reported Locations for Open Space')
legend('3:35pm Open Space','7:35pm Open Space','Google Open Space')
plot_google_map

figure
hold on
plot(lon_house_1,lat_house_1,'.k','MarkerSize',20)
plot(lon_house_2,lat_house_2,'.c','MarkerSize',20)
plot(lon_house_true,lat_house_true,'*r','MarkerSize',20)
hold off
xlabel('Longitude [$^\circ$]')
ylabel('Latitude [$^\circ$]')
title('GPS Application Measurements and Google Maps Reported Locations for Near Building')
legend('3:45pm Building','7:45pm Building','Google Building')
plot_google_map