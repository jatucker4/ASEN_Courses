%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODING CHALLENGE 10 - Integrating Systems of Equations: van der Pol
% 
% The purpose of this challenge is to give you experience formulating
% higher order differential equations using ODE45 with van der Pol's
% equation. van der Pol's equation describes a self-sustaning oscillator,
% essentially a linear spring with a non-linear damper.
%
% Please ZIP and upload your team's script(s) and figures in PNG format to
% Canvas to complete the challenge.
% 
% STUDENT TEAMMATES
% 1.
% 2.
% 3.
% 4.
%
% CHALLENGE AUTHORS
% Allison Anderson, Melinda Zavala
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all
close all
clc

%% Please write your group number and section here:
group = 12;
section = 1;

%% Setting the parameters (given variables)
tspan = [0 20];
mu = 1;
omega0 = 5;
y_dot0 = 0;
y0 = [y_dot0; omega0];

[t, y] = ode45(@(t,y) g_fun(t,y,mu),tspan,y0);
%% Part A

% 1) Plot position vs. velocity, using mu = 1 for 20 seconds. Use the
% initial conditions y(0) = 5 and y_dot(0) = 0.
figure(1)
plot(y(:, 1),y(:,2));
xlabel('position')
ylabel('velocity')
hold on
% 2) On the same figure, plot position vs. velocity with the initial
% conditions y(0) = 2, and y_dot(0) = 0. Use mu = 1 and plot for 20 seconds.
tspan = [0 20];
mu = 1;
omega0 = 2;
y_dot0 = 0;
y0 = [y_dot0; omega0];

[t, y] = ode45(@(t,y) g_fun(t,y,mu),tspan,y0);
plot(y(:, 1),y(:,2));
hold on
% 3) On the same figure, plot position vs. velocity with the initial
% conditions y(0) = 5, and y_dot(0) = 1. Use mu = 1 and plot for 20 seconds.
tspan = [0 20];
mu = 1;
omega0 = 5;
y_dot0 = 1;
y0 = [y_dot0; omega0];

[t, y] = ode45(@(t,y) g_fun(t,y,mu),tspan,y0);
plot(y(:, 1), y(:,2));
% Plot and save figure in PNG format with group number and section in file
% name

saveas(gcf,['mu1_group',num2str(group),'_sec',num2str(section)],'png')

%% Part B
tspan = [0 20];
mu = 1;
omega0 = 5;
y_dot0 = 0;
y0 = [y_dot0; omega0];

[t, y] = ode45(@(t,y) g_fun(t,y,mu),tspan,y0);
% 2) Plot position and velocity from Part A #1, and on the same figure,
% plot position vs. velocity using mu = 2 over *60 seconds*. Use the
% initial conditions y(0) = 5, and y_dot = 0.
figure(2)
plot(y(:, 1),y(:,2));
xlabel('position')
ylabel('velocity')
hold on

tspan = [0 60];
mu = 2;
omega0 = 5;
y_dot0 = 0;
y0 = [y_dot0; omega0];

[t, y] = ode45(@(t,y) g_fun(t,y,mu),tspan,y0);
% Plot and save figure in PNG format with group number and section in file
% name
plot(y(:, 1),y(:,2));
saveas(gcf,['mu2_group',num2str(group),'_sec',num2str(section)],'png')
