%% Housekeeping
clc;
clear all;
close all;
%% Constants
ohm_dot = 1.991063*10^-7;
R = 6378;
T = 1.62*(3600);
r_p = R + linspace(500,600);
r_p = R + 500;
mu = 398600.4415;
J2 = 1.08263*10^-3;
omega_e = (7.2921*10^-5);
%% Start solving for variables
a = (((T/(2*pi))^2) * mu)^(1/3);

e = 1 - (r_p./a);
numerator = (sqrt(mu) * J2 * R^2);
denominator = ((1 - e.^2).^2 * a.^(7/2));

i = acos( -ohm_dot./((3/2) * (numerator./denominator)) );
i = i*(180/pi);

%% Solve for the initial RAAN
lambda_0 = floor((omega_e + ohm_dot)*T * (180/pi));
J_0 = 2458595;
T_0 = (J_0 - 2451545) / 36525;
theta_G0 = 100.4606184 + (36000.77004*T_0) + (0.000387933*(T_0^2)) - ((2.583*(10^-8))*(T_0^3));
multiplier = floor(theta_G0 / 360);
theta_G0 = theta_G0 - multiplier*360;
theta_G = theta_G0 + 360.98564724 * (12/24);
RAAN = theta_G + lambda_0 - 360;

%% Bi-elliptic Hohmann transfer portion
clc; clear; close all;

mu_sun = 132.712e9;
r_Earth = 149.6e6;
r_Mars = 227.9e6;
% Apoapsis radius of the two transfer arcs
r_transfer = 2.7 * (r_Earth);

% First calculate the velocity of the Earth orbiting the Sun
v_A_1 = sqrt(mu_sun/r_Earth);
% Get the semimajor axis of the first transfer ellipse
% Note apoapsis is the r_transfer value and r_periapsis is r_Earth
a_1 = (r_transfer + r_Earth)/2;
v_A_2 = sqrt((2*mu_sun/r_Earth) - mu_sun/a_1);
% Calculate the velocity at the apogee of the first transfer ellipse
v_A_3 = sqrt((2*mu_sun/r_transfer) - mu_sun/a_1);
delta_V_A = abs(v_A_2 - v_A_1);

% First calculate the velocity of the Mars orbiting the Sun
v_B_1 = sqrt(mu_sun/r_Mars);
% Get the semimajor axis of the first transfer ellipse
% Note apoapsis is the r_transfer value and r_periapsis is r_Mars
a_2 = (r_transfer + r_Mars)/2;
v_B_2 = sqrt((2*mu_sun/r_Mars) - mu_sun/a_2);
% Calculate the velocity at the periapsis of the second transfer ellipse
v_B_3 = sqrt((2*mu_sun/r_transfer) - mu_sun/a_2);
delta_V_B = abs(v_B_2 - v_B_1);

% Get the delta v required to go from one transfer ellipse to the other
delta_V_C = abs(v_B_3 - v_A_3);

delta_V_total = delta_V_A + delta_V_B + delta_V_C;

T_Earth = 2*pi*sqrt((r_Earth^3) / mu_sun);
T_Mars = 2*pi*sqrt((r_Mars^3) / mu_sun);
T_transfer_1 = pi*sqrt((a_1^3) / mu_sun);
T_transfer_2 = pi*sqrt((a_2^3) / mu_sun);
t_transfer_total = T_transfer_1 + T_transfer_2;