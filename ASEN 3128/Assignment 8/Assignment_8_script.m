%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment 8 Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
%% Create Constants
k_S = 1:0.01:3;
g = 9.81; % m/s^2
W = (6.366e5) * 4.44822; % N
m = W/g; % kg
rho = 1.2673e-3 * 515.379; % kg/m^3
V = 518 * 0.3048; % m/s
S = 5500 * 0.092903; % m^2;
c_bar = 27.31 * 0.3048; % m
I_x = (1.82e7) * 1.35581795 ;% kg*m^2
I_y = (3.31e7) * 1.35581795 ;% kg*m^2
I_z = (4.97e7) * 1.35581795 ;% kg*m^2
I_zx = (9.7e5) * 1.35581795 ;% kg*m^2
%% Create trim constants
theta_0 = 0;% rad
u_0 = V; % m/s
%% Perform conversions
% Longitudinal
X_u = -4.883e1 * 14.5939029;
X_w = 1.546e3 * 14.5939029;
X_q = 0 * 4.44822;
X_w_dot = 0 * 0.453592;


Z_u = -1.342e3 * 14.5939029;
Z_w = -8.561e3 * 14.5939029;
Z_q = -1.263e5 * 4.44822;
Z_w_dot = 3.104e2 * 0.453592;  


M_u = 8.176e3 * 4.44822;
M_w = -5.627e4 * 4.44822;
M_q = -1.394e7 * 1.3558179483314;
M_w_dot = -4.138e3 * 4.44822;

%% Perform Rotations
% Create DCM (Rotation about Y by -6.8 degrees
zeta = -6.8 * (pi/180);
% Longitudinal Equations
X_u_prime = (X_u * (cos(zeta)^2)) - ((X_w + Z_u)*sin(zeta)*cos(zeta)) + (Z_w*(sin(zeta)^2));
X_w_prime = X_w * (cos(zeta)^2) + (X_u - Z_w)*sin(zeta)*cos(zeta) - Z_u*(sin(zeta)^2);
X_q_prime = X_q * cos(zeta) - Z_q*sin(zeta);
X_w_dot_prime = -Z_w_dot * sin(zeta) * cos(zeta);

Z_u_prime = Z_u * (cos(zeta)^2) - (Z_w - X_u)*sin(zeta)*cos(zeta) - X_w*(sin(zeta)^2);
Z_w_prime = Z_w * (cos(zeta)^2) + (Z_u + X_w)*sin(zeta)*cos(zeta) + X_u*(sin(zeta)^2);
Z_q_prime = Z_q * cos(zeta) + X_q * sin(zeta);
Z_w_dot_prime = Z_w_dot * (cos(zeta)^2);

M_u_prime = M_u * cos(zeta) - M_w * sin(zeta);
M_w_prime = M_w * cos(zeta) + M_u * sin(zeta);
M_q_prime = M_q;
M_w_dot_prime = M_w_dot * cos(zeta);

%% Get the dimensionalized control derivatives
C_x_delta_e = -3.818e-6;
C_z_delta_e = -0.3648;
C_m_delta_e = -1.444;

X_delta_e = C_x_delta_e*0.5*rho*(u_0^2)*S;
Z_delta_e = C_z_delta_e*0.5*rho*(u_0^2)*S;
M_delta_e = C_m_delta_e*0.5*rho*(u_0^2)*S*c_bar;

%% Create the short period approximation matrices and locus plot ie 3b
A_sp = [M_q_prime/I_y, (u_0*M_w_prime)/I_y; 1, 0];

B_sp = [M_delta_e/I_y; 0];

for i = 1:length(k_S)
    k_1_sp = (M_q_prime/M_delta_e)*(1-k_S(i));
    k_2_sp = (u_0*M_w_prime*(1-k_S(i)))/(M_delta_e);
    K_sp = [-k_1_sp, -k_2_sp];
    full_mat_sp = A_sp + B_sp*K_sp;
    evals_sp = eig(full_mat_sp);
    eig_vec_sp(i,:) = evals_sp;
    
    figure(1)
    plot(real(evals_sp),imag(evals_sp),'o')
    hold on
    xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
    ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
    title("$Locus\:Plot\:for\:Short\:Period\:Approx\:Eigenvalues$",'Interpreter',...
        'latex','FontSize',26)
end
%% Construct A matrix
A = [(X_u_prime/m) , (X_w_prime/m) , 0 , -g*cos(theta_0);
     Z_u_prime/(m-Z_w_dot_prime) , Z_w_prime/(m-Z_w_dot_prime), (Z_q_prime + m*u_0)/(m-Z_w_dot_prime), (-m*g*sin(theta_0))/(m-Z_w_dot_prime);
     (1/I_y)*(M_u_prime + (M_w_dot_prime*Z_u_prime)/(m-Z_w_dot_prime)), (1/I_y)*(M_w_prime + (M_w_dot_prime*Z_w_prime)/(m-Z_w_dot_prime)),...
     (1/I_y)*(M_q_prime + (M_w_dot_prime*(Z_q_prime + m*u_0))/(m-Z_w_dot_prime)), -(M_w_dot_prime*m*g*sin(theta_0))/(I_y*(m-Z_w_dot_prime));
     0, 0, 1, 0];
 
evals = eig(A);

%% Construct the B and K matrices as well as the locus plot for full linearized dynamics
B = [X_delta_e/m ; Z_delta_e/(m - Z_w_dot_prime);...
    (M_delta_e/I_y) + (M_w_dot_prime*Z_delta_e)/(I_y*(m - Z_w_dot_prime)); 0];
for i = 1:length(k_S)
    k_1 = (M_q_prime/M_delta_e)*(1-k_S(i));
    k_2 = (u_0*M_w_prime*(1-k_S(i)))/(M_delta_e);
    K = [0, 0, -k_1, -k_2];
    full_mat = A + B*K;
    evals_1 = eig(full_mat);
    eig_vec(i,:) = evals_1;
    
    figure(2)
    plot(real(evals_1),imag(evals_1),'o')
    hold on
    xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
    ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
    title("$Locus\:Plot\:for\:Full\:Longitudinal\:Eigenvalues$",'Interpreter',...
        'latex','FontSize',26)
end

%% Create the first simulation using ks = 1
delta_u_naut = 0;
w_naut = 0;
q_naut = 0;
delta_theta_naut = 0.1;

tspan = [0 1000];

k_1 = (M_q_prime/M_delta_e)*(1-1);
k_2 = (u_0*M_w_prime*(1-1))/(M_delta_e);
K = [0, 0, -k_1, -k_2];
full_mat = A + B*K;
eig_ks1 = eig(full_mat);
y_0 = [delta_u_naut; w_naut; q_naut; delta_theta_naut];
[t,y] = ode45(@(t,y) odefunc(t, y, full_mat), tspan, y_0);
%% Create the plots for the ks = 1 case
figure(3)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta u^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta w^E\:[m/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Full\:Linearized\:Long.\:Simulation\:With\:K_s = 1$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta q\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\theta\:[rad]$",'Interpreter','latex','FontSize',26)
%% Create the first simulation using ks = 2
k_1 = (M_q_prime/M_delta_e)*(1-2);
k_2 = (u_0*M_w_prime*(1-2))/(M_delta_e);
K = [0, 0, -k_1, -k_2];
full_mat = A + B*K;
eig_ks2 = eig(full_mat);
y_0 = [delta_u_naut; w_naut; q_naut; delta_theta_naut];
[t_1,y_1] = ode45(@(t,y) odefunc(t, y, full_mat), tspan, y_0);

%% Create the plots for the ks = 2 case
figure(4)
subplot(2,2,1)
plot(t_1,y_1(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta u^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t_1,y_1(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta w^E\:[m/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Full\:Linearized\:Long.\:Simulation\:With\:K_s = 2$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t_1,y_1(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta q\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t_1,y_1(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\theta\:[rad]$",'Interpreter','latex','FontSize',26)