%% Housekeeping
close all;
clear all;
clc;

%% Create constants
m = 0.068;
g = 9.81;
k = .0024;
r = 0.060;
alpha = (2*10^-6);
zeta = 1*(10^-3);

% Given these relative wind and ang vel vectors
rel_wind = [0;0;0];
ang_vel = [0;0;0.1];

psi = 0;
theta = 0;
phi = 0;

% Then set up the system of equations to solve for the trim thrusts
A = [ 1, 1, -1, -1;...
    -1, 1, 1, -1; ...
    1, -1, 1, -1;...
    1, 1, 1, 1];
b = [(sqrt(2)/r)*alpha*ang_vel(1)*norm(ang_vel);...
    (sqrt(2)/r)*alpha*ang_vel(2)*norm(ang_vel);...
    (1/k)*alpha*ang_vel(3)*norm(ang_vel);...
    m*g*cos(phi)*cos(theta) + zeta*rel_wind(3)*norm(rel_wind)];

trim_thrust_vec = A\b;

f_1 = trim_thrust_vec(1);
f_2 = trim_thrust_vec(2);
f_3 = trim_thrust_vec(3);
f_4 = trim_thrust_vec(4);

% Create the moment of inertia matrix
I_x = 6.8*(10^-5);
I_y = 9.25*(10^-5);
I_z = 1.35*(10^-5);
I_mat = [I_x, 0, 0; 0, I_y, 0; 0, 0, I_z];

% Create the control forces vector(thrust only acts in z)
x_c = 0;
y_c = 0;
z_c = -sum(trim_thrust_vec);
c_forces = [x_c;y_c;z_c];

% Create the control moments vector
L_c = (r/sqrt(2))*(f_1 + f_2 - f_3 - f_4);
M_c = (r/sqrt(2))*(-f_1 + f_2 + f_3 - f_4);
N_c = (k)*(f_1 - f_2 + f_3 - f_4);

c_moments = [L_c; M_c; N_c];

%% Begin the numerical integration sections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State vector:
% [ phi                ]
% [ theta              ]
% [ psi                ] 
% [ p                  ]
% [ q                  ] 
% [ r                  ]
% [ N_inertial_pos     ]
% [ E_inertial_pos     ]
% [ D_inertial_pos     ]
% [ N_inertial_pos_dot ]
% [ E_inertial_pos_dot ]
% [ D_inertial_pos_dot ]
% Describe the initial state and time span for integration:


tspan = [0 5];
% y_0 = [phi, theta, psi, ang_vel(1), ang_vel(2), ang_vel(3), 0, 0, 0, rel_wind(1),...
%     rel_wind(2), rel_wind(3)];
% 
% [t,y] = ode45(@(t,y) odefunc(t,y,c_forces, c_moments,1, 0),...
%     tspan, y_0);

% Uncomment below to use linear model

y_0 = [0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0];

[t,y] = ode45(@(t,y) odefunc(t,y,c_forces, c_moments,0,1),...
    tspan, y_0);
%% Create Plots
% I'm going to plot each state over time
figure(1)
subplot(1,3,1)
plot(t,y(:,3)*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Azimuth\:angle\:[deg]$",'Interpreter','latex','FontSize',26)
xlim([0 5])
% title("$Change\:in\:Azimuth\:Angle\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,2)*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Elevation\:angle\:[deg]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Euler\:Angles\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)
xlim([0 5])

subplot(1,3,3)
plot(t,y(:,1)*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Bank\:angle\:[deg]$",'Interpreter','latex','FontSize',26)
xlim([0 5])
% title("$Change\:in\:Bank\:Angle\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

figure(2)
subplot(1,3,1)
plot(t,y(:,4))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Roll\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Roll\:Rate\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,5))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Pitch\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Angular\:Rates\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(t,y(:,6))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Yaw\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Yaw\:Rate\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

figure(3)
subplot(1,3,1)
plot(t,y(:,7))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$North\:Position\:[m]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:North\:Position\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,8))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$East\:Position\:[m]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Inertial\:Position\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(t,y(:,9))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Upward\:Position\:[m]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Position\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

figure(4)
subplot(1,3,1)
plot(t,y(:,10))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$u\:[m/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:North\:Velocity\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,11))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$v\:[m/s]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Velocity\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(t,y(:,12))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$w\:[m/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Upward\:Velocity\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

% Now a 3D plot of the quad copters position
figure(5)
% y(:,9) = 0;
plot3(y(:,7),y(:,8),y(:,9))
xlabel("$North\:Position\:[m]$",'Interpreter','latex','FontSize',26)
ylabel("$East\:Position\:[m]$",'Interpreter','latex','FontSize',26)
zlabel("$Upward\:Position\:[m]$",'Interpreter','latex','FontSize',26)
title("$Quad\:Copter\:Trajectory$",...
    'Interpreter','latex','FontSize',26)

%% Parse data from the demonstration with assignment gains
S = load('RSdata_White_Assignemnt_Gains.mat');
body_p = S.rt_estim.signals.values(:,10);
body_q = S.rt_estim.signals.values(:,11);
body_r = S.rt_estim.signals.values(:,12);
time_demo = S.rt_estim.time(:);
figure(6)

subplot(1,3,1)
plot(time_demo,body_p)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Roll\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:X\:Position\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(time_demo,body_q)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Pitch\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Angular\:Rates\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(time_demo,body_r)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Yaw\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Z\:Position\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)
%% Plot the experimental drone without assignment gains
S = load('RSdata_White_0903.mat');
body_p = S.rt_estim.signals.values(:,10);
body_q = S.rt_estim.signals.values(:,11);
body_r = S.rt_estim.signals.values(:,12);
time_demo = S.rt_estim.time(:);
figure(7)

subplot(1,3,1)
plot(time_demo,body_p)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Roll\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:X\:Position\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(time_demo,body_q)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Pitch\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Angular\:Rates\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(time_demo,body_r)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Yaw\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Z\:Position\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)