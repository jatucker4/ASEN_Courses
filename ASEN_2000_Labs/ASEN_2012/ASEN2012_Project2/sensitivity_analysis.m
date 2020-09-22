%% Housekeeping
clear;
clc;
close all;
global m_bottle_temp C_drag_temp theta_temp Vol_water_initial
%% First vary the mass of the water
water_vol_temp = linspace(0.0003,0.001,100);
for i = 1:length(water_vol_temp)

    %% Create variables
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
    m_bottle_temp = 0.112;
    C_drag_temp = 0.33;
    P_initial = 275790 + P_atm;
    Vol_water_initial = water_vol_temp(i);
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
    %% Modify for ISP model
    % First process ISP data
    % [isp,delta_v] = isp_process(6.1,40,977.1,1.6520,'Group_10_08AM_statictest2.txt');
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
    range_vec(i) = range;
    height_vec(i) = height;
    %% Plot Verification Case
%     figure(1)
%     plot3(x(1:index),y(1:index),z(1:index));
%     xlabel('$Horizontal\:Distance\:[m]$',"Interpreter","latex","FontSize",...
%         15);
%     ylabel('$Cross\:Distance\:[m]$',"Interpreter","latex","FontSize",...
%         15);
%     zlabel('$Vertical\:Distance\:[m]$',"Interpreter","latex","FontSize",...
%         15);
%     title('$3D\:Thermodynamic\:Predicted\:Trajectory$',...
%         "Interpreter","latex","FontSize",15);
%     legend("Predicted Trajectory");

end