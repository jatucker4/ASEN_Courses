%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver script for the second part of Lab 1 in ASEN 3200 and
% is used to solve Kepler's equation iteratively using Newton-Raphson.
% Elliptical orbits only.
%
% Created by: Johnathan Tucker
%
% Date Created: 9/3/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all;
close all;
clc;
%% Create the necessary constants
% Generally this problems are begin with a known time and period
t = 13.5 * 60; % [s]
T = 22680; % [s]
% And if position is to be calculated then eccentricity and semi-major axis
% will be needed
e = 0.52972;
a = 3997.57; % [km]
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

%% Now calculate the position if necessary
r = a*(1-e*cos(E)); % [km]
% Check to see if this position is within the observation range
r_moon = 1737.1; % [km]
r_obs = r - r_moon; % [km]