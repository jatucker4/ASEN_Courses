%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is for Homework 5 problem 3.12
%
% Created by: Johnathan Tucker
%
% Created on: 2/15/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
%% Create necessary constants for part a
theta = linspace(pi,0.5);
x = ((pi-theta)./sin(theta)).*cos(theta);
y = (pi-theta);
x(1) = -1;
%% Create plot for part a
figure(1)
plot(x,y)
title("Resulting Semi-Infinite Body")
xlabel("x-axis")
ylabel("y-axis")
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%% Create cp function for part b
cp = -((1./((x.^2)+(y.^2))) + (2.*x./((x.^2)+(y.^2))));
figure(2)
plot(x,cp)
title("Coefficient of Pressure Plot")
xlabel("x-axis")
ylabel("Cp")
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';