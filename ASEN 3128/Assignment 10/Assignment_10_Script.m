%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment 10 Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
%% Create Constants
zeta = -6.8; % degrees
g = 9.81; % m/s^2
W = (6.366e5) * 4.44822; % N
m = W/g; % kg
rho = 1.2673e-3 * 515.379; % kg/m^3
V = 518 * 0.3048; % m/s
S = 5500 * 0.092903; % m^2;
b = 195.68 * 0.3048; % m
c_bar = 27.31 * 0.3048; % m
I_x = (1.82e7) * 1.35581795 ;% kg*m^2
I_y = (3.31e7) * 1.35581795 ;% kg*m^2
I_z = (4.97e7) * 1.35581795 ;% kg*m^2
I_zx = (9.7e5) * 1.35581795 ;% kg*m^2

% Rotate inertia values to the stability frame
I_x_r = I_x*(cosd(zeta)^2) + I_z*(sind(zeta)^2) + I_zx*sind(2*zeta);% kg*m^2
I_y_r = I_y;% kg*m^2
I_z_r = I_x*(sind(zeta)^2) + I_z*(cosd(zeta)^2) - I_zx*sind(2*zeta);% kg*m^2
I_zx_r = -1*(0.5*(I_x - I_z)*sind(2*zeta) + I_zx*(sind(zeta)^2 - cosd(zeta)^2));% kg*m^2

I_x_p = (I_x_r*I_z_r - I_zx_r^2)/I_z_r;% kg*m^2
I_z_p = (I_x_r*I_z_r - I_zx_r^2)/I_x_r;% kg*m^2
I_zx_p = I_zx_r/(I_x_r*I_z_r - I_zx_r^2);% kg*m^2


% Create trim constants
u_0 = V;
theta_0 = 0;

%% Create non-dimensional stability derivative constants
C_y_B = -0.8771;
C_y_p = 0;
C_y_r = 0;

C_l_B = -0.2797;
C_l_p = -0.3295;
C_l_r = 0.304;

C_n_B = 0.1964;
C_n_p = -0.04073;
C_n_r = -0.2737;
%% Calculate dimensional stability derivatives in SI units
Y_v = 0.5*rho*u_0*S*C_y_B;
Y_p = 0.25*rho*u_0*b*S*C_y_p;
Y_r = 0.25*rho*u_0*b*S*C_y_r;

L_v = 0.5*rho*u_0*b*S*C_l_B;
L_p = 0.25*rho*u_0*(b^2)*S*C_l_p;
L_r = 0.25*rho*u_0*(b^2)*S*C_l_r;

N_v = 0.5*rho*u_0*b*S*C_n_B;
N_p = 0.25*rho*u_0*(b^2)*S*C_n_p;
N_r = 0.25*rho*u_0*(b^2)*S*C_n_r;
%% Create New Table for Unit Converted Stability Derivatives
row_names ={'_';'u[m/s]'; 'p[rad/s]'; 'r[rad/s]'};
Y_vec = {'Y[N]';Y_v; Y_p; Y_r};
L_vec = {'L[Nm]'; L_v; L_p; L_r};
N_vec = {'N[Nm]'; N_v; N_p; N_r};
T = table(Y_vec, L_vec, N_vec);
T.Properties.RowNames = row_names;
table2latex(T,'non_stab_table');
%% Convert to the Stability Frame
Y_v_p = Y_v;
Y_p_p = Y_p*cosd(zeta) - Y_r*sind(zeta);
Y_r_p = Y_r*cosd(zeta) + Y_p*sind(zeta);

L_v_p = L_v*cosd(zeta) - N_v*sind(zeta);
L_p_p = L_p*(cosd(zeta)^2) - (L_r + N_p)*sind(zeta)*cosd(zeta) + N_r*(sind(zeta)^2);
L_r_p = L_r*(cosd(zeta)^2) - (N_r - L_p)*sind(zeta)*cosd(zeta) - N_p*(sind(zeta)^2);

N_v_p = N_v*cosd(zeta) + L_v*sind(zeta);
N_p_p = N_p*(cosd(zeta)^2) - (N_r - L_p)*sind(zeta)*cosd(zeta) - L_r*(sind(zeta)^2);
N_r_p = N_r*(cosd(zeta)^2) + (L_r + N_p)*sind(zeta)*cosd(zeta) + L_p*(sind(zeta)^2);
%% Create New Table for Unit Converted Stability Derivatives- Stability frame
row_names ={'_';'u[m/s]'; 'p[rad/s]'; 'r[rad/s]'};
Y_vec = {'Y[N]';Y_v_p; Y_p_p; Y_r_p};
L_vec = {'L[Nm]'; L_v_p; L_p_p; L_r_p};
N_vec = {'N[Nm]'; N_v_p; N_p_p; N_r_p};
T = table(Y_vec, L_vec, N_vec);
T.Properties.RowNames = row_names;
table2latex(T,'stab_table');
%% Build the A matrix
A = [Y_v_p/m, Y_p_p/m, (Y_r_p/m) - u_0, g*cos(theta_0);...
    (L_v_p/I_x_p) + I_zx_p*N_v_p, (L_p_p/I_x_p) + I_zx_p*N_p_p, (L_r_p/I_x_p) + I_zx_p*N_r_p, 0;...
    (I_zx_p*L_v_p + (N_v_p/I_z_p)), (I_zx_p*L_p_p + (N_p_p/I_z_p)), (I_zx_p*L_r_p + (N_r_p/I_z_p)), 0;...
    0, 1, tan(theta_0), 0];
%% Get the eigenvalues and eigenvectors for the A matrix
[evecs, evals] = eig(A);
%% Pull out the different modes and plot them
DR_mode = [evals(1,1),evals(2,2)];
roll_mode = evals(3,3);
spiral_mode = evals(4,4);
figure(50)
plot(real(DR_mode),imag(DR_mode),'o')
hold on
plot(roll_mode,0,'o')
hold on
plot(spiral_mode,0,'o')
    xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
    ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
    title("$Locus\:Plot\:for\:Lateral\:Eigenvalues$",'Interpreter',...
        'latex','FontSize',26)
    legend("Dutch Roll Mode","Roll Mode","Spiral Mode")
    
%% Get the time constants for each
% Dutch Roll Mode
omega_d_dutch = imag(evals(1,1));
nat_freq_dutch = sqrt(omega_d_dutch^2 + real(evals(1,1))^2);
damping_ratio_dutch = -real(evals(1,1))/nat_freq_dutch;
time_const_dutch = 1/(damping_ratio_dutch*nat_freq_dutch);
dutch_period = 1/damping_ratio_dutch;
dutch_period = (2*pi)/omega_d_dutch;
% Rolling Convergence
time_const_rolling = -1/evals(3,3);

% Spiral Mode
time_const_spiral = -1/evals(4,4);
%% Create New Table for the modal information
row_names ={'Mode:';'Dutch Roll'; 'Roll'; 'Spiral'};
Y_vec = {'\tau [s]';time_const_dutch; time_const_rolling; time_const_spiral};
L_vec = {'\omega_n [rad/s]'; nat_freq_dutch; '-'; '-'};
N_vec = {'\zeta'; damping_ratio_dutch; '-'; '-'};
T = table(Y_vec, L_vec, N_vec);
T.Properties.RowNames = row_names;
table2latex(T,'modal_info');
%% Create the dutch roll mode approximation
% For reference the dutch roll calculated from the A matrix is:
% -0.0149 pm 0.9136i
fancy_Y_v = A(1,1);
fancy_N_v = A(3,1);
fancy_N_r = A(3,3);

dutch_approx_1 = 0.5*(fancy_Y_v + fancy_N_r)+ 0.5*sqrt((fancy_Y_v + fancy_N_r)^2 - 4*(fancy_Y_v*fancy_N_r + u_0*fancy_N_v));
dutch_approx_2 = 0.5*(fancy_Y_v + fancy_N_r)- 0.5*sqrt((fancy_Y_v + fancy_N_r)^2 - 4*(fancy_Y_v*fancy_N_r + u_0*fancy_N_v));

% Get the damping ratio, natural freq, and period of approx for comparison
omega_d_dutch_approx = imag(dutch_approx_1);
nat_freq_dutch_approx = sqrt(omega_d_dutch_approx^2 + real(dutch_approx_1)^2);
damping_ratio_dutch_approx = -real(dutch_approx_1)/nat_freq_dutch_approx;
time_const_dutch_approx = 1/(damping_ratio_dutch_approx*nat_freq_dutch_approx);
dutch_period_approx = 1/damping_ratio_dutch_approx;
dutch_period_approx = (2*pi)/omega_d_dutch_approx;
% Get percent difference between approx and actual for damping ratio,
% natural freq, and period
perc_diff_damping_ratio = 100*(abs(damping_ratio_dutch - damping_ratio_dutch_approx)/damping_ratio_dutch_approx);
perc_diff_nat_freq = 100*(abs(nat_freq_dutch-nat_freq_dutch_approx)/nat_freq_dutch_approx);
perc_diff_period = 100*(abs(dutch_period-dutch_period_approx)/dutch_period_approx);

fprintf("The percent difference between the dutch roll damping ratio and approximation is: %f\n ",perc_diff_damping_ratio);
fprintf("The percent difference between the dutch roll natural frequency and approximation is: %f\n ",perc_diff_nat_freq);
fprintf("The percent difference between the dutch roll period and approximation is: %f\n ",perc_diff_period);

%% Start the simulation by veryfying the trim condition
delta_v_0 = 0;
delta_p_0 = 0;
delta_r_0 = 0;
delta_phi_0 = 0;

tspan = [0 50];
y_0 = [delta_v_0; delta_p_0; delta_r_0; delta_phi_0];
[t,y] = ode45(@(t,y) ode_func(t, y, A), tspan, y_0);

% Create plots
figure(1)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta v^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta p^E\:[rad/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Verifying\:Trim\:Condition$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta r\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\phi\:[rad]$",'Interpreter','latex','FontSize',26)

%% Deviation of delta_v_0
delta_v_0 = 10;
delta_p_0 = 0;
delta_r_0 = 0;
delta_phi_0 = 0;

tspan = [0 100];
y_0 = [delta_v_0; delta_p_0; delta_r_0; delta_phi_0];
[t,y] = ode45(@(t,y) ode_func(t, y, A), tspan, y_0);

% Create plots
figure(2)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta v^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta p^E\:[rad/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Initial\:Deviation\:of\:\Delta v^E = 10[m/s]$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta r\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\phi\:[rad]$",'Interpreter','latex','FontSize',26)

%% Response to Deviation of delta_p_0
delta_v_0 = 0;
delta_p_0 = 0.15;
delta_r_0 = 0;
delta_phi_0 = 0;

tspan = [0 100];
y_0 = [delta_v_0; delta_p_0; delta_r_0; delta_phi_0];
[t,y] = ode45(@(t,y) ode_func(t, y, A), tspan, y_0);

% Create plots
figure(4)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta v^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta p^E\:[rad/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Initial\:Deviation\:of\:\Delta p=0.15[rad/s]$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta r\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\phi\:[rad]$",'Interpreter','latex','FontSize',26)

%% Response to given y(o) = [-1.8563,-0.4185,0.0311,0.6148]T
tspan = [0 100];
y_0 = [-1.8563,-0.4185,0.0311,0.6148];
[t,y] = ode45(@(t,y) ode_func(t, y, A), tspan, y_0);

% Create plots
figure(5)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta v^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta p^E\:[rad/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Response\:to\:Initial\:Condition:[-1.8563,-0.4185,0.0311,0.6148]^T$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta r\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\phi\:[rad]$",'Interpreter','latex','FontSize',26)

%% Response to given y(o) = [2.9477,-0.0063,0.0758,1.2431]T
tspan = [0 100];
y_0 = [2.9477,-0.0063,0.0758,1.2431];
[t,y] = ode45(@(t,y) ode_func(t, y, A), tspan, y_0);

% Create plots
figure(6)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta v^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta p^E\:[rad/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Response\:to\:Initial\:Condition:[2.9477,-0.0063,0.0758,1.2431]^T$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta r\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\phi\:[rad]$",'Interpreter','latex','FontSize',26)
%% Project the given initial vectors onto the eigenspace to determine which modes are activated
y_init_1 = [-1.8563,-0.4185,0.0311,0.6148]';
y_init_2 = [2.9477,-0.0063,0.0758,1.2431]';

proj_1 = evecs*y_init_1;
proj_2 = evecs*y_init_2;