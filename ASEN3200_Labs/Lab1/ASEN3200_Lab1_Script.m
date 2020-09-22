%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver script for the first part of Lab 1 in ASEN 3200
%
% Created by: Johnathan Tucker
%
% Date Created: 9/3/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all;
close all;
clc;
%% Define constants
a = 242200; % [km]
e = 0.587;
x_0 = 100029; % [km]
y_0 = 0; % [km]
z_0 = 0; % [km]
xdot_0 = 0; % [km/s]
ydot_0 = 2.47655; % [km/s]
zdot_0 = 0.436682; % [km/s]
state_init = [x_0, y_0, z_0, xdot_0, ydot_0, zdot_0];
tspan = [0 1186238.974];
%% Plug into ode45
[t,y] = ode45(@(t,y) odefunc(t,y), tspan, state_init);
%% Create plots
figure(1)
plot3(y(:,1),y(:,2),y(:,3))
hold on
scatter3(0,0,0)
hold on
scatter3(x_0,y_0,z_0)
hold on
quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[100029;0;0],[0;1;0],[0;0;1]);
title("Satellite Orbit Around the Earth")
legend("Satellite Orbit", "Earth", "Satellite initial position",...
    "Eccentricity Vector")
xlabel("x-axis (km)")
ylabel("y-axis (km)")
zlabel("z-axis (km)")
%% Get the magnitude of the eccentricity and angular momentum 
t_vec = linspace(0,1186238.974);
e_vec = ones(1,length(t_vec)).*e;
% Recall: r_p = a(1-e) and h = r*v_perp but the satellite is at perigee so
% all velocity is in the perp direction
r_p = a*(1-e);
h = norm([xdot_0,ydot_0,zdot_0])*r_p;
h_vec = ones(1,length(t_vec)).*h;
%% Now plot them over the span
figure(2)
plot(t_vec,e_vec)
title("Magnitude Eccentricity of Over One Orbital Period")
legend("Eccentricity Magnitude")
xlabel("Time(s)")
ylabel("Unitless")

figure(3)
plot(t_vec,h_vec)
title("Magnitude of Anglular Momentum Over One Orbital Period")
legend("Angular Momentum Magnitude")
xlabel("Time(s)")
ylabel("h (km^3/s^2")