%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment 7 Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
%% Create Constants
g = 9.81;
W = (6.366e5) * 4.44822;
m = W/g;
V = 518 * 0.3048;
I_x = (1.82e7) * 1.35581795 ;
I_y = (3.31e7) * 1.35581795 ;
I_z = (4.97e7) * 1.35581795 ;
I_zx = (9.7e5) * 1.35581795 ;
%% Create trim constants
theta_0 = 0;
u_0 = V;
%% Perform conversions
% Longitudinal
X_u = -4.883e1 * 14.5939029;
X_w = 1.546e3 * 14.5939029;
X_q = 0 * 4.44822;
X_w_dot = 0 * 0.453592;
X_delta_e = 3.994e4 * 4.44822;

Z_u = -1.342e3 * 14.5939029;
Z_w = -8.561e3 * 14.5939029;
Z_q = -1.263e5 * 4.44822;
Z_w_dot = 3.104e2 * 0.453592;  
Z_delta_e = -3.341e5 * 4.44822;

M_u = 8.176e3 * 4.44822;
M_w = -5.627e4 * 4.44822;
M_q = -1.394e7 * 1.3558179483314;
M_w_dot = -4.138e3 * 4.44822;
M_delta_e = -3.608e7 * 1.3558179483314;

%% Create New Table for Unit Converted Stability Derivatives
row_names ={'_';'u[m/s]'; 'w[m/s]'; 'q[rad/s]'; '\Dot{w}[ft/s^2]'; '\delta_e[rad]'};
X_vec = {'X[N]';X_u; X_w; X_q; X_w_dot; X_delta_e};
Z_vec = {'Z[N]'; Z_u; Z_w; Z_q; Z_w_dot; Z_delta_e};
M_vec = {'M[Nm]'; M_u; M_w; M_q; M_w_dot; M_delta_e};
T = table(X_vec, Z_vec, M_vec);
T.Properties.RowNames = row_names;
table2latex(T);
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
%% Create New Table for Rotation Converted Stability Derivatives
row_names ={'_';'u[m/s]'; 'w[m/s]'; 'q[rad/s]'; '\Dot{w}[ft/s^2]'; '\delta_e[rad]'};
X_vec = {'X[N]';X_u_prime; X_w_prime; X_q_prime; X_w_dot_prime; X_delta_e};
Z_vec = {'Z[N]'; Z_u_prime; Z_w_prime; Z_q_prime; Z_w_dot_prime; Z_delta_e};
M_vec = {'M[Nm]'; M_u_prime; M_w_prime; M_q_prime; M_w_dot_prime; M_delta_e};
T = table(X_vec, Z_vec, M_vec);
T.Properties.RowNames = row_names;
table2latex(T,'table2.tex');
%% Construct A matrix
A = [(X_u_prime/m) , (X_w_prime/m) , 0 , -g*cos(theta_0);
     Z_u_prime/(m-Z_w_dot_prime) , Z_w_prime/(m-Z_w_dot_prime), (Z_q_prime + m*u_0)/(m-Z_w_dot_prime), (-m*g*sin(theta_0))/(m-Z_w_dot_prime);
     (1/I_y)*(M_u_prime + (M_w_dot_prime*Z_u_prime)/(m-Z_w_dot_prime)), (1/I_y)*(M_w_prime + (M_w_dot_prime*Z_w_prime)/(m-Z_w_dot_prime)),...
     (1/I_y)*(M_q_prime + (M_w_dot_prime*(Z_q_prime + m*u_0))/(m-Z_w_dot_prime)), -(M_w_dot_prime*m*g*sin(theta_0))/(I_y*(m-Z_w_dot_prime));
     0, 0, 1, 0];
 
evals = eig(A);
%% Corresponding natural frequencies and damping ratios
% First short mode evals
real_short = real(evals(1));
omega_d_short = imag(evals(1));

short_omega_n = sqrt(omega_d_short^2 + real_short^2);
short_zeta = real_short/-short_omega_n;

% Now phugoid mode evals
real_phu = real(evals(3));
omega_d_phu = imag(evals(3));

phu_omega_n = sqrt(omega_d_phu^2 + real_phu^2);
phu_zeta = real_phu/-phu_omega_n;
%% Short period mode approximation
% Real part first:
lambda_real = M_q_prime/(2*I_y);
% Imaginary part:
lambda_imaginary = (1/(2*I_y))*sqrt(M_q_prime^2 + 4*I_y*u_0*M_w_prime);

%% Lanchester Approximation period
lan_approx = pi*sqrt(2)*u_0/g;
eval_phu_period = (2*pi)/imag(evals(3));

%% Deviation of delta_u_naut = 10
delta_u_naut = 10;
w_naut = 0;
q_naut = 0;
delta_theta_naut = 0;

tspan = [0 1000];
y_0 = [delta_u_naut; w_naut; q_naut; delta_theta_naut];
[t,y] = ode45(@(t,y) linearized_long(t, y, A), tspan, y_0);

%% Create plots
figure(1)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta u^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta w^E\:[m/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Initial\:Deviation\:of\:\Delta u^E = 10[m/s]$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta q\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\theta\:[rad]$",'Interpreter','latex','FontSize',26)

%% Deviation of w_naut = 10
delta_u_naut = 0;
w_naut = 10;
q_naut = 0;
delta_theta_naut = 0;

tspan = [0 1000];
y_0 = [delta_u_naut; w_naut; q_naut; delta_theta_naut];
[t,y] = ode45(@(t,y) linearized_long(t, y, A), tspan, y_0);

%% Create plots
figure(2)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta u^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta w^E\:[m/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Initial\:Deviation\:of\:\Delta w^E =10[m/s]$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta q\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\theta\:[rad]$",'Interpreter','latex','FontSize',26)

%% Deviation of q_naut = 0.1
delta_u_naut = 0;
w_naut = 0;
q_naut = 0.1;
delta_theta_naut = 0;

tspan = [0 1000];
y_0 = [delta_u_naut; w_naut; q_naut; delta_theta_naut];
[t,y] = ode45(@(t,y) linearized_long(t, y, A), tspan, y_0);

%% Create plots
figure(3)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta u^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta w^E\:[m/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Initial\:Deviation\:of\:\Delta q=0.1[rad/s]$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta q\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\theta\:[rad]$",'Interpreter','latex','FontSize',26)

%% Deviation of delta_theta_naut = 0.1
delta_u_naut = 0;
w_naut = 0;
q_naut = 0;
delta_theta_naut = 0.1;

tspan = [0 1000];
y_0 = [delta_u_naut; w_naut; q_naut; delta_theta_naut];
[t,y] = ode45(@(t,y) linearized_long(t, y, A), tspan, y_0);

%% Create plots
figure(4)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta u^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta w^E\:[m/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Initial\:Deviation\:of\:\Delta\theta=0.1[rad]$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta q\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\theta\:[rad]$",'Interpreter','latex','FontSize',26)
%% Trim state plots
delta_u_naut = 0;
w_naut = 0;
q_naut = 0;
delta_theta_naut = 0;

tspan = [0 1000];
y_0 = [delta_u_naut; w_naut; q_naut; delta_theta_naut];
[t,y] = ode45(@(t,y) linearized_long(t, y, A), tspan, y_0);


figure(5)
subplot(2,2,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta u^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta w^E\:[m/s]$",'Interpreter','latex','FontSize',26)
sgtitle("$Trim\:State\:$",...
    'Interpreter','latex','FontSize',26)

subplot(2,2,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta q\:[rad/s]$",'Interpreter','latex','FontSize',26)

subplot(2,2,4)
plot(t,y(:,4))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta\theta\:[rad]$",'Interpreter','latex','FontSize',26)