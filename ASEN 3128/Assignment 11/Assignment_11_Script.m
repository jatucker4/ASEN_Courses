%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment 11 Script
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
row_names ={'_';'$\delta a$'; '$\delta r$[rad/s]'};
Y_vec = {'Y[N]';Y_delta_a; Y_delta_r;};
L_vec = {'L[Nm]';L_delta_a; L_delta_r;};
N_vec = {'N[Nm]';N_delta_a; N_delta_r;};
T = table(Y_vec, L_vec, N_vec);
T.Properties.RowNames = row_names;
table2latex(T,'cont_table');
%% Create the 4x2 B matrix
B = [Y_delta_a/m , Y_delta_r/m ; ...
    (L_delta_a/I_x_p) + I_zx_p*N_delta_a ,(L_delta_r/I_x_p) + I_zx_p*N_delta_r;...
    I_zx_p*L_delta_a + (N_delta_a/I_z_p), I_zx_p*L_delta_r + (N_delta_r/I_z_p);...
    0, 0];
%% Create augmented B and A matrices
A_aug = [A, zeros(4,2); 0, 0, sec(theta_0), 0, 0, 0;1, 0, 0, 0, u_0*cos(theta_0),0];
B_aug = [B; zeros(2,2)];

%% Begin question 3
% Part a pos feedback
k_vec = 0:0.01:10;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, 0, k_vec(i), 0, 0; 0, 0, 0, 0, 0, 0];
    full_mat = A_aug + B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(1)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
hold on
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3a$",'Interpreter',...
        'latex','FontSize',26)


% Part b neg feedback
k_vec = 0:0.01:10;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, k_vec(i), 0, 0, 0, 0; 0, 0, 0, 0, 0, 0];
    full_mat = A_aug - B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(2)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3b$",'Interpreter',...
        'latex','FontSize',26)

% Part c neg feedback
k_vec = 0:0.01:10;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, k_vec(i), 0, 0, 0; 0, 0, 0, 0, 0, 0];
    full_mat = A_aug - B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(3)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3c$",'Interpreter',...
        'latex','FontSize',26)
    
% Part d pos feedback
k_vec = 0:0.01:20;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, 0, 0, k_vec(i), 0; 0, 0, 0, 0, 0, 0];
    full_mat = A_aug + B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(4)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3d$",'Interpreter',...
        'latex','FontSize',26)

% Part e neg feedback
k_vec = 0:0.0001:0.1;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, 0, 0, 0, 0; k_vec(i), 0, 0, 0, 0, 0];
    full_mat = A_aug - B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(5)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3e$",'Interpreter',...
        'latex','FontSize',26)

% Part f neg feedback
k_vec = 0:0.01:2;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, 0, 0, 0, 0; 0, k_vec(i), 0, 0, 0, 0];
    full_mat = A_aug - B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(6)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3f$",'Interpreter',...
        'latex','FontSize',26)
    
% Part g pos feedback
k_vec = 0:0.01:5;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, 0, 0, 0, 0; 0, 0, k_vec(i), 0, 0, 0];
    full_mat = A_aug + B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(7)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3g$",'Interpreter',...
        'latex','FontSize',26)

% Part g neg feedback
k_vec = 0:0.01:5;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, 0, 0, 0, 0; 0, 0, k_vec(i), 0, 0, 0];
    full_mat = A_aug - B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(8)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3g\:Negative\:Feedback$",'Interpreter',...
        'latex','FontSize',26)


% Part h neg feedback
k_vec = 0:0.01:5;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, 0, 0, 0, 0; 0, 0, 0, k_vec(i), 0, 0];
    full_mat = A_aug - B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(9)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3h$",'Interpreter',...
        'latex','FontSize',26)
    
% Part h pos feedback
k_vec = 0:0.01:5;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, 0, 0, 0, 0; 0, 0, 0, k_vec(i), 0, 0];
    full_mat = A_aug + B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(10)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3h\:Positive\:Feedback$",'Interpreter',...
        'latex','FontSize',26)

% Part i pos feedback
k_vec = 0:0.01:5;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, k_vec(i), 0];
    full_mat = A_aug + B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(11)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3i$",'Interpreter',...
        'latex','FontSize',26)
    
% Part i neg feedback
k_vec = 0:0.01:5;
for i = 1:length(k_vec)
    % Create U vector
    K_mat = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, k_vec(i), 0];
    full_mat = A_aug - B_aug*K_mat;
    evals = eig(full_mat);
%     real_eval_vec(:,i) = real(evals);
%     imag_eval_vec(:,i) = imag(evals);
    figure(12)
    if i == 1
        p1 = plot(real(evals),imag(evals),'r*');
    elseif i == length(k_vec)
        p3 = plot(real(evals),imag(evals),'b*');
    else
        plot(real(evals),imag(evals),'k.')
    end
    hold on
end
legend([p1 p3],{'Minimum Gains','Maximum Gains'})
xlabel("$Real\:Axis$",'Interpreter','latex','FontSize',26);
ylabel("$Imaginary\:Axis$",'Interpreter','latex','FontSize',26);
title("$Locus\:Plot\:for\:Lateral\:Eigenvalues\:Case\:3i\:Negative\:Feedback$",'Interpreter',...
        'latex','FontSize',26)