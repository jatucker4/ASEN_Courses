clear;
close all;
clc;
N = 1000;
global m_bottle_temp theta_temp C_drag_temp P_initial Vol_water_initial wind_surface wind_aloft out_mat
m_bottle = 0.112 + .001*randn(1,N);
theta = (45 + 1*randn(1,N)) * (pi/180);
C_drag = 0.33 + 0.005*randn(1,N);
vol_water = 0.0005 + (5*10^(-7))*randn(1,N);
press_temp = 1*randn(1,N);
wind_surf_temp = 2 + 1*randn(1,N);
out_mat = sts_dataprocess('Group_10_08AM_statictest2.txt');
for i = 1:N
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
    m_bottle_temp = m_bottle(i);
    C_drag_temp = C_drag(i);
%     P_initial = 275790 + P_atm;
    Vol_water_initial = vol_water(i);
    T_initial = 300;
    v_0 = 0;
    theta_temp = theta(i);
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
    [wind_surface,wind_aloft] = getWind([wind_surf_temp(i),...
        deg2rad(-22.5),0],[0,deg2rad(22.5),0]);
    
    Vol_air_initial = Vol_bottle - Vol_water_initial;
    m_water_init = Vol_water_initial * rho_water;
    m_air_initial = (P_initial*Vol_air_initial)/(R*T_initial);
    m_initial_total = m_water_init + m_air_initial + m_bottle_temp;
    %% Modify for ISP model
%     out_mat = sts_dataprocess('Group_10_08AM_statictest2.txt');
    isp = out_mat{1,7};
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
    crossrange = y_ode(index);
    x_mc(i) = range;
    y_mc(i) = crossrange;
    x_mat{i,:} = x(1:index);
    y_mat{i,:} = y_ode(1:index);
    z_mat{i,:} = z(1:index);
end

x = x_mc;
y = y_mc;
figure(1);
plot(x,y, "k.","markersize",6)
axis equal;
grid on;
xlabel('$Horizontal\:Distance\:[m]$',"Interpreter","latex","FontSize",...
    15);
ylabel('$Cross\:Range\:Distance\:[m]$',"Interpreter","latex","FontSize",...
    15);
title('$Error\:Ellipse\:With\:Projected\:Monte\:Carlo\:Landing\:Spots$',...
    "Interpreter","latex","FontSize",15);

P = cov(x,y);
mean_x = mean(x);
mean_y = mean(y);

n = 100;
p = 0:pi/n:2*pi;

[eigvec,eigval] = eig(P);
xy_vect = [cos(p'), sin(p')] * sqrt(eigval) * eigvec';
x_vect = xy_vect(:,1);
y_vect = xy_vect(:,2);
hold on
plot(1*x_vect+mean_x,1*y_vect+mean_y,"b");
hold on
plot(2*x_vect+mean_x,2*y_vect+mean_y,"g");
hold on
plot(3*x_vect+mean_x,3*y_vect+mean_y,"r");
hold on

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
Vol_water_initial = 0.0005;
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
% [wind_aloft,wind_surface] = getWind([2,deg2rad(-22.5),0],[0,deg2rad(22.5),0]);
%% Modify for ISP model

%     out_mat = sts_dataprocess('Group_10_08AM_statictest2.txt');
    isp = out_mat{1,7};
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
    crossrange = y(index);
% plot(range,height,"*","markersize",9,"b*")
plot(range,crossrange,'color', [.6 .2 1],'linestyle','none',...
    'marker','x','markersize',15,'LineWidth',2);
legend("Simulated Landing Spots","1\sigma","2\sigma","3\sigma",...
    "Predicted Landing Spot");
for i = 1:length(z_mat)
    x_temp = x_mat{i};
    y_temp = y_mat{i};
    z_temp = z_mat{i};
    figure(2)
    hold on
    plot3(x_temp,y_temp,z_temp)
end
xlabel('$Horizontal\:Distance\:[m]$',"Interpreter","latex","FontSize",...
    15);
ylabel('$Cross\:Distance\:[m]$',"Interpreter","latex","FontSize",...
    15);
zlabel('$Vertical\:Distance\:[m]$',"Interpreter","latex","FontSize",...
    15);
title('$Projected\:Monte\:Carlo\:Trajectories$',...
    "Interpreter","latex","FontSize",15);