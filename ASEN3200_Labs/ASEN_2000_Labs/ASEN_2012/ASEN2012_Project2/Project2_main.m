%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main script for the second project in ASEN 2012 
%
% Created by: Johnathan Tucker
% Author ID: 108602275
%
% Purpose: The purpose of this script is to generate a thrust and
% trajectory profile for a water bottle rocket
%
% Inputs: N/A
%
% Outputs: Trajectory and thrust profile plots
%
% Assumptions: 
%   - Static stability
%   - No wind
%   - Forces due to rail are ignored
%
% Date Created:  11/23/2018
% Date Modified: 11/24/2018
%                11/26/2018
%                11/28/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear;
clc;
close all;
global m_bottle_temp C_drag_temp theta_temp P_initial wind_aloft wind_surface Vol_water_initial
%% Create variables
isp = 1;
g = 9.81;
C_discharge = 0.8;
rho_air_amb = 0.961;
rho_air_amb = 1.184;
Vol_bottle = 0.002;
P_atm = 83054.2464;
gamma = 1.4;
rho_water = 1000;
d_throat = .021;
d_bottle = .105;
R = 287;
m_bottle_temp = 0.117;
C_drag_temp = 0.33;
P_initial = 275790 + P_atm;
Vol_water_initial = 0.00962;
T_initial = 300;
v_0 = 0;
theta_temp = 45 * (pi/180);
x_0 = 0;
z_0 = 0.25;
tspan = [0 20];
rail_length = 0.5;
vel_x_init = 0;
vel_y_init = 0;
vel_z_init = 0;
y_0 = 0;
A_throat = pi*(d_throat/2)^2;
A_bottle = pi*(d_bottle/2)^2;

Vol_air_initial = Vol_bottle - Vol_water_initial;

m_water_init = Vol_water_initial * rho_water;
m_air_initial = (P_initial*Vol_air_initial)/(R*T_initial);
m_initial_total = m_water_init + m_air_initial + m_bottle_temp;
%% Modify for ISP model
% First process ISP data
if isp
    out_mat = sts_dataprocess('Group_10_08AM_statictest2.txt');
    isp = out_mat{1,7};
    [wind_aloft,wind_surface] = getWind([1,deg2rad(22.5),0],[3,deg2rad(22.5),0]);
    delta_V = isp*9.81*log(m_initial_total/m_bottle_temp);
    vel_Vec = [cos(theta_temp); 0 ; sin(theta_temp)] .* delta_V;
    vel_x_init = vel_Vec(1);
    vel_y_init = vel_Vec(2);
    vel_z_init = vel_Vec(3);
    P_initial = P_atm;
    initial_conditions = [x_0, y_0,z_0,vel_x_init, vel_y_init,...
        vel_z_init,0,0,Vol_bottle];
    [t, y] = ode45('Project_2_diffeqs',tspan,initial_conditions);
    x = y(:,1);
    z = y(:,3);
    y_ode = y(:,2);
    % Find Range and height of Rocket
    z_0 = 0;
    % Look for when rocket hits ground
    hits_ground = find(x > 30, 1, 'first');
    diff_vec = abs(z(hits_ground:end) - z_0);
    % Find index of the closest value
    [~, index] = min(diff_vec);
    % Adjust the index
    index = index + length(z(1:hits_ground - 1)); 
    height = max(z);
    range = x(index);
    
    %% Plot Verification Case
    figure(1)
    plot3(x(1:index),y(1:index),z(1:index));
    xlabel('$Horizontal\:Distance\:[m]$',"Interpreter","latex","FontSize",...
        15);
    ylabel('$Cross\:Distance\:[m]$',"Interpreter","latex","FontSize",...
        15);
    zlabel('$Vertical\:Distance\:[m]$',"Interpreter","latex","FontSize",...
        15);
    title('$3D\:ISP\:Predicted\:Trajectory$',...
        "Interpreter","latex","FontSize",15);
    legend("Predicted Trajectory");
    % Get thrust vector from function
    % thrust = thrust_profile(y);
    % 
    % figure(2)
    % plot(t,thrust)
    % xlim([0,0.45])
    % ylim([0,200])
    % xlabel('Time (s)')
    % ylabel('Thrust (N)')
    % title('Rocket Thrust versus Time')
    fprintf('Max height the rocket achieved is %0.3g (m)\n',height);
    fprintf('Max distance the rocket achieved is %0.3g (m)\n',range);
else
%% Create initial conditions
initial_conditions = [x_0, y_0,z_0,vel_x_init, vel_y_init,...
    vel_z_init,m_air_initial,m_water_init,Vol_air_initial];

%% Perform ODE calculations
[t, y] = ode45('Project_2_diffeqs',tspan,initial_conditions);

x = y(:,1);
z = y(:,3);
y = y(:,2);
% Find Range and height of Rocket
z_0 = 0;
% Look for when rocket hits ground
hits_ground = find(x > 30, 1, 'first');
diff_vec = abs(z(hits_ground:end) - z_0);
% Find index of the closest value
[~, index] = min(diff_vec);
% Adjust the index
index = index + length(z(1:hits_ground - 1)); 
height = max(z);
range = x(index);

%% Plot Verification Case
figure(1)
plot3(x(1:index),y(1:index),z(1:index));
xlabel('$Horizontal\:Distance\:[m]$',"Interpreter","latex","FontSize",...
    15);
ylabel('$Cross\:Distance\:[m]$',"Interpreter","latex","FontSize",...
    15);
zlabel('$Vertical\:Distance\:[m]$',"Interpreter","latex","FontSize",...
    15);
title('$3D\:Thermodynamic\:Predicted\:Trajectory$',...
    "Interpreter","latex","FontSize",15);
legend("Predicted Trajectory");
% Get thrust vector from function
% thrust = thrust_profile(y);
% 
% figure(2)
% plot(t,thrust)
% xlim([0,0.45])
% ylim([0,200])
% xlabel('Time (s)')
% ylabel('Thrust (N)')
% title('Rocket Thrust versus Time')
fprintf('Max height the rocket achieved is %0.3g (m)\n',height);
fprintf('Max distance the rocket achieved is %0.3g (m)\n',range);
end
%% Ellipse finding stuff
