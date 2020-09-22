%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver script for the sixth question in homework 5
%
% Created by: Johnathan Tucker
%
% Date Created: 10/7/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all;
close all;
clc;
%% Create constants
mu = 132712000000;
r = 1.433*(10^9);
e = 0;
T = 9.356*(10^8);
t = 2.75846*(10^8);
% Calculate the mean anomaly
M_e = (t/T) * 2*pi; % [rad]

%% Begin implementing Newton-Raphson
% This error was chosen based off numerical analysis standards. This is how
% many digits MATLAB is precise to and therefore it's our tolerance.
tolerance = 1*(10^(-16));

% As per the book on page 153 choose an initial guess for E
if M_e < pi
E = M_e + e/2;
else
E = M_e - e/2;
end
% This initial value of delta is necessary for the while loop to start
delta = 1;
% Every iteration check if the difference is within the tolerance
while abs(delta) > tolerance
    delta = (E - e*sin(E) - M_e)/(1 - e*cos(E));
    E = E - delta;
end

%% Get perifocal position for Saturn
r_s = [r*cos(E), r*sin(E)];

% Now get perifocal position for Earth
r_2 = 149.6*(10^6);
