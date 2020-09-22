%% Housekeeping
clc;
clear all;
close all;
%% Create constants
k1 = 0.001276;
k2 = 0.00232;
k3_vec = linspace(-0.5,0.5,30000);
Ix = 5.8*10^-5;
g = -9.81;
%% Create the system matrix
for i = 1:length(k3_vec)
    k3 = k3_vec(i);
    A = [-k1/Ix, -k2/Ix, -(k2*k3)/Ix ; 1, 0, 0; 0, g, 0];
    evals = eig(A);
    flag1 = -1/evals(1);
    flag2 = -1/evals(2);
    flag3 = -1/evals(3);
    if ((flag1>0 && flag1 <= 1.25 && isreal(flag1)) && (flag2>0 && flag2<= 1.25 && isreal(flag2)) && (flag3>0 && flag3<= 1.25 && isreal(flag3)))
        eval_vec(:,i) = evals;
        figure(1)
    end    
end
 % I get that K3 = -0.046984899496650
 Iy = 7.2*10^-5;
for i = 1:length(k3_vec)
    k3 = k3_vec(i);
    A = [-k1/Iy, -k2/Iy, -(k2*k3)/Iy ; 1, 0, 0; 0, g, 0];
    evals = eig(A);
    flag1 = -1/evals(1);
    flag2 = -1/evals(2);
    flag3 = -1/evals(3);
    if ((flag1>0 && flag1 <= 1.25 && isreal(flag1)) && (flag2>0 && flag2<= 1.25 && isreal(flag2)) && (flag3>0 && flag3<= 1.25 && isreal(flag3)))
        eval_vec_2(:,i) = evals;
        figure(1) 
    end    
end
%% Create the Locus plots for each:
% Lateral first
k3_vec = linspace(-2,2);
for i = 1:length(k3_vec)
    k3 = k3_vec(i);
    A = [-k1/Ix, -k2/Ix, -(k2*k3)/Ix ; 1, 0, 0; 0, g, 0];
    evals = eig(A);
    figure(1)
    plot(real(evals),imag(evals),'o');
    hold on
    xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
    ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
    title("$Locus\:Plot\:for\:Lateral\:Eigenvalues$",'Interpreter',...
        'latex','FontSize',26)
end
for i = 1:length(k3_vec)
    k3 = k3_vec(i);
    A = [-k1/Iy, -k2/Iy, -(k2*k3)/Iy ; 1, 0, 0; 0, g, 0];
    evals = eig(A);
    figure(2)
    plot(real(evals),imag(evals),'o');
    hold on
    xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
    ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
    title("$Locus\:Plot\:for\:Longitudinal\:Eigenvalues$",'Interpreter',...
        'latex','FontSize',26)
end
%% Setup for question 2
m = 0.068;
g = 9.81;
k = .0024;
r = 0.060;
alpha = (2*10^-6);
zeta = 1*(10^-3);

% Given these relative wind and ang vel vectors
rel_wind = [0;0;0];
ang_vel = [0;0;0];

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
% Solve for the trim thrusts
trim_thrust_vec = A\b;

f_1 = trim_thrust_vec(1);% trim thrust [N]
f_2 = trim_thrust_vec(2);% trim thrust [N]
f_3 = trim_thrust_vec(3);% trim thrust [N]
f_4 = trim_thrust_vec(4);% trim thrust [N]

% Create the moment of inertia matrix
I_x = 5.8*(10^-5);% moment of inertia [kg m^2]
I_y = 7.2*(10^-5);% moment of inertia [kg m^2]
I_z = 1.0*(10^-5);% moment of inertia [kg m^2]
I_mat = [I_x, 0, 0; 0, I_y, 0; 0, 0, I_z];% moment of inertia [kg m^2]

% Create the control forces vector(thrust only acts in z)
x_c = 0; % control force x [N]
y_c = 0;% control force y [N]
z_c = -sum(trim_thrust_vec);% control force x [N]
c_forces = [x_c;y_c;z_c];

% Create the control moments vector
L_c = (r/sqrt(2))*(f_1 + f_2 - f_3 - f_4);% control moment x [Nm]
M_c = (r/sqrt(2))*(-f_1 + f_2 + f_3 - f_4);% control moment y [Nm]
N_c = (k)*(f_1 - f_2 + f_3 - f_4);% control moment z [Nm]

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
tspan = [0 2];
dev_u_ref = 0;
dev_v_ref = 0.5;
y_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

[t,y] = ode45(@(t,y) odefunclateral(t,y,c_forces, c_moments,0,1,dev_v_ref,...
    dev_u_ref),tspan, y_0);


%% Create Plots for the lateral case
% I'm going to plot each state over time
figure(3)
subplot(1,3,1)
plot(t,y(:,3)*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Azimuth\:angle\:[deg]$",'Interpreter','latex','FontSize',26)

% title("$Change\:in\:Azimuth\:Angle\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,2)*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Elevation\:angle\:[deg]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Euler\:Angles\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)


subplot(1,3,3)
plot(t,y(:,1)*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Bank\:angle\:[deg]$",'Interpreter','latex','FontSize',26)

% title("$Change\:in\:Bank\:Angle\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

figure(4)
subplot(1,3,1)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Roll\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Roll\:Rate\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,5))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Pitch\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Angular\:Rates\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(t,y(:,6))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Yaw\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Yaw\:Rate\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

figure(5)
subplot(1,3,1)
plot(t,y(:,7))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$North\:Position\:[m]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:North\:Position\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,8))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$East\:Position\:[m]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Inertial\:Position\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(t,y(:,9))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Upward\:Position\:[m]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Position\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

figure(6)
subplot(1,3,1)
plot(t,y(:,10))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$u\:[m/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:North\:Velocity\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,11))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$v\:[m/s]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Velocity\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(t,y(:,12))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$w\:[m/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Upward\:Velocity\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

% Now a 3D plot of the quad copters position
figure(7)
% y(:,9) = 0;
plot3(y(:,7),y(:,8),y(:,9))
xlabel("$North\:Position\:[m]$",'Interpreter','latex','FontSize',26)
ylabel("$East\:Position\:[m]$",'Interpreter','latex','FontSize',26)
zlabel("$Upward\:Position\:[m]$",'Interpreter','latex','FontSize',26)
title("$Quad\:Copter\:Trajectory$",...
    'Interpreter','latex','FontSize',26)
%% Solve for the Longitudinal reference velocity
dev_v_ref = 0;
tolerance = 1*10^-3;
y_0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

[t,y] = ode45(@(t,y) odefunclongitudinal(t,y,c_forces, c_moments,0,1,dev_v_ref,...
    0.5),tspan, y_0);

%% Create Plots for the longitudinal case
% I'm going to plot each state over time
figure(8)
subplot(1,3,1)
plot(t,y(:,3)*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Azimuth\:angle\:[deg]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Azimuth\:Angle\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,2)*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Elevation\:angle\:[deg]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Euler\:Angles\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(t,y(:,1)*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Bank\:angle\:[deg]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Bank\:Angle\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

figure(9)
subplot(1,3,1)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Roll\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Roll\:Rate\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,5))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Pitch\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Angular\:Rates\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(t,y(:,6))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Yaw\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Yaw\:Rate\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

figure(10)
subplot(1,3,1)
plot(t,y(:,7))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$North\:Position\:[m]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:North\:Position\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,8))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$East\:Position\:[m]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Inertial\:Position\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(t,y(:,9))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Upward\:Position\:[m]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Position\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

figure(11)
subplot(1,3,1)
plot(t,y(:,10))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$u\:[m/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:North\:Velocity\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(t,y(:,11))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$v\:[m/s]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Velocity\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(t,y(:,12))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$w\:[m/s]$",'Interpreter','latex','FontSize',26)
% title("$Change\:in\:Upward\:Velocity\:Over\:Time$",...
%     'Interpreter','latex','FontSize',26)

% Now a 3D plot of the quad copters position
figure(12)
% y(:,9) = 0;
plot3(y(:,7),y(:,8),y(:,9))
xlabel("$North\:Position\:[m]$",'Interpreter','latex','FontSize',26)
ylabel("$East\:Position\:[m]$",'Interpreter','latex','FontSize',26)
zlabel("$Upward\:Position\:[m]$",'Interpreter','latex','FontSize',26)
title("$Quad\:Copter\:Trajectory$",'Interpreter','latex','FontSize',26)
%% NOTES
% The eigenvalues are close together which causes the sluggish response. No
% one eigenvalue is dominating. Reference the locus as to what behavior the
% evals you chose obey ie underdamped, overdamped, critically damped(last
% one) 
%% Plot the Experimental results
% Load the data struct
S = load('RSdata_one_1013.mat');
% Pull the necessary data from the struct
body_p = smooth(S.rt_estimatedStates.signals.values(:,10));
body_q = smooth(S.rt_estimatedStates.signals.values(:,11));
body_r = smooth(S.rt_estimatedStates.signals.values(:,12));
body_u = S.rt_estimatedStates.signals.values(:,7);
body_v = S.rt_estimatedStates.signals.values(:,8);
body_w = S.rt_estimatedStates.signals.values(:,9);
body_phi = S.rt_estimatedStates.signals.values(:,6);
body_theta = S.rt_estimatedStates.signals.values(:,5);
body_psi = S.rt_estimatedStates.signals.values(:,4);
time_demo = S.rt_estimatedStates.time(:);

% Begin plotting
figure(13)
subplot(1,3,1)
plot(time_demo,body_p)
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Roll\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(time_demo,body_q)
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Pitch\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Angular\:Rates\:Over\:Time\:Mambo\:Quad\:Copter$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(time_demo,body_r)
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Yaw\:Rate\:[rad/s]$",'Interpreter','latex','FontSize',26)


figure(14)
subplot(1,3,1)
plot(time_demo,body_u)
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$u\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(1,3,2)
plot(time_demo,body_v)
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$v\:[m/s]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Velocity\:Over\:Time\:Mambo\:Quad\:Copter$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(time_demo,body_w)
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$w\:[m/s]$",'Interpreter','latex','FontSize',26)


figure(15)
subplot(1,3,1)
plot(time_demo,body_psi*(180/pi))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Azimuth\:angle\:[deg]$",'Interpreter','latex','FontSize',26)


subplot(1,3,2)
plot(time_demo,body_theta*(180/pi))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Elevation\:angle\:[deg]$",'Interpreter','latex','FontSize',26)
title("$Change\:in\:Euler\:Angles\:Over\:Time\:Mambo\:Quad\:Copter$",...
    'Interpreter','latex','FontSize',26)

subplot(1,3,3)
plot(time_demo,body_psi*(180/pi))
xlim([0 5])
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Bank\:angle\:[deg]$",'Interpreter','latex','FontSize',26)