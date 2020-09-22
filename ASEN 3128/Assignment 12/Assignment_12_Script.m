%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment 12 Script
%
%
% This script creates all plots and calls the ode functions
%
% Created by: Johnathan Tucker
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
% row_names ={'_';'u[m/s]'; 'p[rad/s]'; 'r[rad/s]'};
% Y_vec = {'Y[N]';Y_v; Y_p; Y_r};
% L_vec = {'L[Nm]'; L_v; L_p; L_r};
% N_vec = {'N[Nm]'; N_v; N_p; N_r};
% T = table(Y_vec, L_vec, N_vec);
% T.Properties.RowNames = row_names;
% table2latex(T,'non_stab_table');
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
% row_names ={'_';'u[m/s]'; 'p[rad/s]'; 'r[rad/s]'};
% Y_vec = {'Y[N]';Y_v_p; Y_p_p; Y_r_p};
% L_vec = {'L[Nm]'; L_v_p; L_p_p; L_r_p};
% N_vec = {'N[Nm]'; N_v_p; N_p_p; N_r_p};
% T = table(Y_vec, L_vec, N_vec);
% T.Properties.RowNames = row_names;
% table2latex(T,'stab_table');
%% Build the A matrix
A = [Y_v_p/m, Y_p_p/m, (Y_r_p/m) - u_0, g*cos(theta_0);...
    (L_v_p/I_x_p) + I_zx_p*N_v_p, (L_p_p/I_x_p) + I_zx_p*N_p_p, (L_r_p/I_x_p) + I_zx_p*N_r_p, 0;...
    (I_zx_p*L_v_p + (N_v_p/I_z_p)), (I_zx_p*L_p_p + (N_p_p/I_z_p)), (I_zx_p*L_r_p + (N_r_p/I_z_p)), 0;...
    0, 1, tan(theta_0), 0];
evals_non_aug = eig(A);
%% Create fancy variables
fancy_Y_v = A(1,1);
fancy_Y_r = A(1,3);
fancy_L_v = A(2,1);
fancy_L_p = A(2,2);
fancy_L_r = A(2,3);
fancy_N_v = A(3,1);
fancy_N_p = A(3,2);
fancy_N_r = A(3,3);

%% Create Spiral mode approximation
E = g*((fancy_L_v*fancy_N_r - fancy_L_r*fancy_N_v)*cos(theta_0) + ...
    (fancy_L_p*fancy_N_v - fancy_L_v*fancy_N_p)*sin(theta_0));
D = -g*(fancy_L_v*cos(theta_0) + fancy_N_v*sin(theta_0)) +...
    fancy_Y_v*(fancy_L_r*fancy_N_p - fancy_L_p*fancy_N_r) + ...
    fancy_Y_r*(fancy_L_p*fancy_N_v - fancy_L_v*fancy_N_p);

spiral_approx = -E/D;

%% Create rolling mode approximation
rolling_approx = fancy_L_p;
% Create a locus plot comparing the approximations to the actual values
figure(15)
plot(rolling_approx,0,'b*')
hold on
plot(spiral_approx,0,'r*')
hold on
plot(real(evals_non_aug(4,:)),imag(evals_non_aug(4,:)),'c*')
hold on
plot(real(evals_non_aug(3,:)),imag(evals_non_aug(3,:)),'g*')
legend('Roll Mode Approx.','Spiral Mode Approx', 'Spiral Mode','Roll Mode')
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:Comparing\:Lateral\:Roll\:and\:Spiral\:Mode\:Approximations$",'Interpreter',...
        'latex','FontSize',26)
%% Finding the lateral control derivatives
% Non-dimensional
C_y_delta_a = 0;
C_y_delta_r = .1146;

C_L_delta_a = -1.368e-2;
C_L_delta_r = 6.976e-3;

C_N_delta_a = -1.973e-4;
C_N_delta_r = -0.1257;

% Dimensional
Y_delta_a = 0.5*rho*(u_0^2)*S*C_y_delta_a;
Y_delta_r = 0.5*rho*(u_0^2)*S*C_y_delta_r;

L_delta_a = 0.5*rho*(u_0^2)*S*b*C_L_delta_a;
L_delta_r = 0.5*rho*(u_0^2)*S*b*C_L_delta_r;
 
N_delta_a = 0.5*rho*(u_0^2)*S*b*C_N_delta_a;
N_delta_r = 0.5*rho*(u_0^2)*S*b*C_N_delta_r;
%% Create New Table for Control Derivatives
% row_names ={'_';'$\delta a$'; '$\delta r$[rad/s]'};
% Y_vec = {'Y[N]';Y_delta_a; Y_delta_r;};
% L_vec = {'L[Nm]';L_delta_a; L_delta_r;};
% N_vec = {'N[Nm]';N_delta_a; N_delta_r;};
% T = table(Y_vec, L_vec, N_vec);
% T.Properties.RowNames = row_names;
% table2latex(T,'cont_table');
%% Create the 4x2 B matrix
B = [Y_delta_a/m , Y_delta_r/m ; ...
    (L_delta_a/I_x_p) + I_zx_p*N_delta_a ,(L_delta_r/I_x_p) + I_zx_p*N_delta_r;...
    I_zx_p*L_delta_a + (N_delta_a/I_z_p), I_zx_p*L_delta_r + (N_delta_r/I_z_p);...
    0, 0];
%% Create augmented B and A matrices
A_aug = [A, zeros(4,2); 0, 0, sec(theta_0), 0, 0, 0;1, 0, 0, 0, u_0*cos(theta_0),0];
B_aug = [B; zeros(2,2)];

%% Get natural time constants for the spiral and DR modes
evals = eig(A_aug);
time_const_spiral = -1/evals(5);
time_const_DR = -1/real(evals(3));

omega_d_dutch = imag(evals(3,1));
nat_freq_dutch = sqrt(omega_d_dutch^2 + real(evals(3,1))^2);
damping_ratio_dutch = -real(evals(3,1))/nat_freq_dutch;
%% I'm going to create the requirement thresholds for the control design
% Dutch roll damping ratio
zeta_DR = 0.35; % Must be greater than or equal to this
slope = tan(acos(zeta_DR));
% Dutch roll real component
real_DR = -0.025; % Must be greater than this

% Spiral real component
real_SP = -0.04; % Must be greater than this

% Initial condition vector
y_0 = [10, -.14, .05, 0, 0, 0]';

% Based on assignment 11 I'd like to try the K del_R and K del_p schemes
% The gain values in each of these were:
K_del_r = 0:0.01:3;
K_del_p = 0:.01:10; %0:0.01:1.4;

for i = 1:length(K_del_r)
    % Create U vector
    K_mat = [0, 0, 0, 0, 0, 0; 0, 0, K_del_r(i), 0, 0, 0];
    full_mat = A_aug + B_aug*K_mat;
    evals = eig(full_mat);
    figure(1)
    if length(find(imag(evals)))> 2
        p3 = plot(real(evals),imag(evals),'b*');
        hold on
        break
    end
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
        hold on
    elseif i == length(K_del_r)
        p3 = plot(real(evals),imag(evals),'b*');
        hold on
    else
        plot(real(evals),imag(evals),'k.');
        hold on
    end
end
p4 = plot(real(evals),-slope*real(evals));
hold on
p5 = plot(real(evals),slope*real(evals));
hold on
p6 = xline(-0.04);
hold on
p7 = xline(-0.025);
legend([p1 p3 p4 p5 p6 p7],{'Minimum Gains','Maximum Gains',...
    'Cone Boundary','Cone Boundary','Spiral Real Boundary',...
    'Dutch Roll Real Boundary'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:of\:Boundaries\:and\:Chosen\:Control\:Law$",'Interpreter',...
        'latex','FontSize',26)
% From the plot I can choose k_del_r to be 1.43/1.24
%% Check the behavior over time using ode45
% Create constants
K_mat = [0, 0, 0, 0, 0, 0; 0, 0, 1.43, 0, 0, 0];
tspan = [0 100];
[t,y] = ode45(@(t,y) ode_func(t, y, A_aug,B_aug,K_mat) , tspan, y_0);

% Plot the delta v over time and delta psi over time to check overshoot
figure(6)
subplot(2,1,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta v^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,1,2)
plot(t,y(:,5).*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta \psi\:[deg]$",'Interpreter','latex','FontSize',26)
sgtitle("$Closed\:Loop\:System\:Response$",...
    'Interpreter','latex','FontSize',26)
fprintf("Minimum overshoot of delta v is: %f\n",min(y(:,1)))
fprintf("Maximum overshoot of delta v is: %f\n",max(y(:,1)))
fprintf("Minimum overshoot of delta psi is: %f\n",min(y(:,5)).*(180/pi))
fprintf("Maximum overshoot of delta psi is: %f\n",max(y(:,5)).*(180/pi))


% Calculate max deflection angles
for i = 1:length(y)
    temp_2(:,i) = K_mat*y(i,:)';
end
delta_a = temp_2(1,:);
delta_r = temp_2(2,:);
% Write the max deflection angle to the command line in degrees
fprintf("Maximum deflection of delta a is: %f\n",max(delta_a).*(180/pi))
fprintf("Maximum deflection of delta r is: %f\n",max(delta_r).*(180/pi))

% Plot the aileron and rudder deflection over time
figure(16)
subplot(2,1,1)
plot(t,delta_a.*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta \delta a\:[deg]$",'Interpreter','latex','FontSize',26)
sgtitle("$Aileron\:and\:Rudder\:Deflection\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

subplot(2,1,2)
plot(t,delta_r.*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta \delta r\:[deg]$",'Interpreter','latex','FontSize',26)

%% Perform other checks using A_aug, B_aug and the gains on natural behavior
check_mat = A_aug + B_aug*K_mat;
check_evals = eig(check_mat);
figure(10)
cp1 = plot(real(check_evals),imag(check_evals),'o');
hold on
cp2 = plot(real(check_evals),slope*real(check_evals));
hold on
cp3 = plot(real(check_evals),-slope*real(check_evals));
hold on
cp4 = xline(-0.04);
hold on
cp5 = xline(-0.025);
legend([cp1,cp2,cp3,cp4,cp5],...
    {'Eigenvalues','Cone Boundary', 'Cone Boundary','Spiral Real Boundary',...
    'Dutch Roll Real Boundary'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Closed\:Loop\:Eigenvalue\:Target\:Region\:With\:Chosen\:Gain\:Eigenvalues$",'Interpreter',...
        'latex','FontSize',26)

figure(11)
cp21 = plot(real(check_evals),slope*real(check_evals));
hold on
cp31 = plot(real(check_evals),-slope*real(check_evals));
hold on
cp41 = xline(-0.04);
hold on
cp51 = xline(-0.025);
legend([cp21,cp31,cp41,cp51],...
    {'Cone Boundary', 'Cone Boundary','Spiral Real Boundary',...
    'Dutch Roll Real Boundary'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Closed\:Loop\:Eigenvalue\:Target\:Region$",'Interpreter',...
        'latex','FontSize',26)
% Calculate the time constants for each mode as well as damping ratio of DR
% to verify
omega_d_dutch_check = imag(check_evals(3,1));
nat_freq_dutch_check  = sqrt(omega_d_dutch_check^2 + real(check_evals(3,1))^2);
damping_ratio_dutch_check  = -real(check_evals(3,1))/nat_freq_dutch_check;

check_DR_time_const = -1/real(check_evals(3,:));
check_spiral_time_const = -1/check_evals(5);
check_roll_time_const = -1/check_evals(6);

fprintf("Dutch Roll damping ratio is: %f\n",damping_ratio_dutch_check)
fprintf("Dutch Roll time constant is: %f\n",check_DR_time_const)
fprintf("Roll time constant is: %f\n",check_roll_time_const)
fprintf("Spiral time constant is: %f\n",check_spiral_time_const)

%% Create plots for every closed loop state response
figure(2)
subplot(2,3,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta v^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,3,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta p\:[rad/s]$",'Interpreter','latex','FontSize',26)


subplot(2,3,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta r\:[rad/s]$",'Interpreter','latex','FontSize',26)


subplot(2,3,4)
plot(t,y(:,4).*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta \phi\:[deg]$",'Interpreter','latex','FontSize',26)


subplot(2,3,5)
plot(t,y(:,5).*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta \psi\:[rad]$",'Interpreter','latex','FontSize',26)


subplot(2,3,6)
plot(t,y(:,6))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta y^E\:[m]$",'Interpreter','latex','FontSize',26)
sgtitle("$Closed\:Loop\:System\:Response$",...
    'Interpreter','latex','FontSize',26)

%% Now use ode45 to get the natural system response
[t,y] = ode45(@(t,y) ode_func_2(t, y, A_aug) , tspan, y_0);
figure(3)
subplot(2,3,1)
plot(t,y(:,1))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta v^E\:[m/s]$",'Interpreter','latex','FontSize',26)

subplot(2,3,2)
plot(t,y(:,2))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta p\:[rad/s]$",'Interpreter','latex','FontSize',26)


subplot(2,3,3)
plot(t,y(:,3))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta r\:[rad/s]$",'Interpreter','latex','FontSize',26)


subplot(2,3,4)
plot(t,y(:,4).*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta \phi\:[deg]$",'Interpreter','latex','FontSize',26)

subplot(2,3,5)
plot(t,y(:,5).*(180/pi))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta \psi\:[rad]$",'Interpreter','latex','FontSize',26)

subplot(2,3,6)
plot(t,y(:,6))
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$\Delta y^E\:[m]$",'Interpreter','latex','FontSize',26)
sgtitle("$Open\:Loop\:System\:Response$",...
    'Interpreter','latex','FontSize',26)


