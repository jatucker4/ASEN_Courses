%% Housekeeping
clc;
clear all;
close all;
%% Create Mach and Cp vectors
M_vec = linspace(0,1,1000);
Cp = -0.41./(sqrt(1-M_vec.^2));
Cp_cr = (2./(1.4.*M_vec.^2)).*((((1+(.4./2).*M_vec.^2)./(1+(.4./2))).^(1.4/.4)) - 1);

%% Find the intersection
[~,index] = min(abs(Cp-Cp_cr));
Mach_cr = M_vec(index)
