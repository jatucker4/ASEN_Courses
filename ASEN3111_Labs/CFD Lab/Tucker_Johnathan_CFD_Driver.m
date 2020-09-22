%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3111 - CFD
% 
% Created By: Johnathan Tucker
%
% Collaborators: 
%
% The purpose of the script is to act as a driver that will execute the
% functions necesary to solve questions two and three of CA4. Note that the
% solution to question 1 is the PLLT function itself
%
% Created Date: 4/9/2020
%
% Change Log: 
%           - 4/7/2020 Code the PLLT function
%           - 4/8/2020 Code up questions two and three
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
tic
%Global formatting commands to imporve graphing looks:
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex'); 
%% First put the experimental CFD data into vectors
alpha_CFD = [-5:1:12 , 14, 16];
cl_CFD = [-.32438, -.21503, -.10081, .010503, .12155, .24163, .34336,...
    .45256, .56037, .66625, .76942, .86923, .96386, 1.0441, 1.0743,...
    1.0807, 1.0379, 1.034, 1.0156, 0.97946];
cd_CFD = [.044251, .033783, .028627, .025864, .024643, .025099, .025635,...
    .02766, .030677, .034855, .040403, .04759, .057108, .070132, .090921,...
    .11193, .13254, .15645, .20959, .25668];

%% Next put the MH 32 data into vectors
alpha_2d = -5:1:15;
CL_2d = [-0.2446, -0.1465, -.0401, .0568, .1717, .2737, .4058, .5143,...
    .6167, .7194, .8201, .9193, 1.0129, 1.1027, 1.1844, 1.2533, 1.2865,...
    1.2763, 1.2329, 1.1635, 1.0951];
CD_2d = [0.0140, .0091, .0073, .0059, .0049, .0043, .0045, .0050, .0057,...
    .0066, .0078, .0092, .0112, .0134, .0165, .0201, .0252, .0332,...
    .0475, .072, .1052];
%% Create necessary constants and feed them into the PLLT function to get CL and CD
b = 3.22; % [m]
% Fit a line to the MH 32 data to get the cross sectional lift slope
% ensuring to only fit to the linear region.
p = polyfit(alpha_2d(1:16),CL_2d(1:16),1);
a0_t = p(1)*(180/pi);
a0_r = p(1)*(180/pi);
c_t = .08;
c_r = .23;
aero_t = interp1(CL_2d,alpha_2d,0);
aero_r = aero_t;
N = 500;
for i = 1:length(alpha_2d)
    [~,c_L_PLLT(i),c_Di_PLLT(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,alpha_2d(i),alpha_2d(i),N);
end

figure
subplot(1,3,1)
hold on
plot(alpha_CFD,cl_CFD)
plot(alpha_2d,c_L_PLLT)
hold off
grid on
xlabel("$\alpha$ [deg]","FontSize",14)
ylabel("$C_l$","FontSize",14)
legend("Synthetic","PLLT")
title("Tempest/Twister Synthetic $C_l$ vs $\alpha$","FontSize",14)

subplot(1,3,2)
plot(alpha_CFD,cd_CFD)
grid on
xlabel("$\alpha$ [deg]","FontSize",14)
ylabel("$C_D$","FontSize",14)
title("Tempest/Twister Synthetic $C_D$ vs $\alpha$","FontSize",14)

subplot(1,3,3)
plot(alpha_CFD,cl_CFD./cd_CFD)
grid on
xlabel("$\alpha$ [deg]","FontSize",14)
ylabel("$\frac{L}{D}$","FontSize",14)
title("Tempest/Twister Synthetic $\frac{L}{D}$ vs $\alpha$","FontSize",14)

%% Start the second question
% Construct the linear fit for the Synthetic data
% First fit a line to get the slope
p_synth = polyfit(alpha_CFD(1:13),cl_CFD(1:13),1);
% Create a line using the slope and y intercept
cl_synth_line = p_synth(1).*(-5:1:7) + p_synth(2);
% Plot the results
figure
hold on
plot(alpha_CFD,cl_CFD)
plot(-5:1:7,cl_synth_line)
hold off
grid on
xlabel("$\alpha$ [deg]","FontSize",14)
ylabel("$C_l$","FontSize",14)
legend("Synthetic","Linear Regression")
title("Tempest/Twister Synthetic $C_l$ vs $\alpha$ with Linear Regression","FontSize",14)

% Estimate for Oswalds
e = 1.78*(1 - 0.045*(16.5^.68)) - 0.64;
K = 1/(16.5*pi*e);
% Interpolate CD0 for when L = 0 //0.025099
CD0_CFD = interp1(cl_CFD,cd_CFD,0);
[CDmin_CFD,index] = min(cd_CFD);

% First calculate the drag polar using 1.1 and 1.2
CD_1_1 = CDmin_CFD + K.*(cl_CFD(1:13) -cl_CFD(index)).^2;
CD_1_2 = CD0_CFD + (cl_CFD(1:13).^2)./(pi*e*16.5);
% Plot the CFD Drag polar for the linear portion of the lift curve
figure
hold on
plot(cl_CFD(1:13),cd_CFD(1:13))
plot(cl_CFD(1:13),CD_1_1)
plot(cl_CFD(1:13),CD_1_2)
grid on
xlabel("$C_L$","FontSize",14)
ylabel("$C_D$","FontSize",14)
title("Tempest/Twister Drag Polar Using Synthetic Data, EQN 1.1, and EQN 1.2","FontSize",14)
legend("Synthetic","EQN 1.1","EQN 1.2")

%% Begin calculations for the fourth question
rho_inf = 1.9869e-3 * 515.379; %[kg/m^3]
V_inf = 18; % [m/s]
W = 9.81*(6.4); %[N]
S = 0.63; %[m^2]
cl_cruise = W/(0.5*rho_inf*S*V_inf^2);
% Interpolate the CD values that correspond to this CL value for each
% calculation method
cd_interp_synth = interp1(cl_CFD,cd_CFD,cl_cruise);
cd_interp_1_1 = interp1(cl_CFD(1:13),CD_1_1,cl_cruise);
cd_interp_1_2 = interp1(cl_CFD(1:13),CD_1_2,cl_cruise);
% Print the CD values closest to this CL value
fprintf("The operating point on the Synthetic drag polar is CL = %f and CD = %f\n",cl_cruise,cd_interp_synth)
fprintf("The operating point on the EQN 1.1 drag polar is CL = %f and CD = %f\n",cl_cruise,cd_interp_1_1)
fprintf("The operating point on the EQN 1.2 drag polar is CL = %f and CD = %f\n\n",cl_cruise,cd_interp_1_2)

% Show on the drag polar plot where the operting points are
figure
hold on
plot(cl_CFD(1:13),cd_CFD(1:13))
plot(cl_CFD(1:13),CD_1_1)
plot(cl_CFD(1:13),CD_1_2)
scatter(cl_cruise,cd_interp_synth)
scatter(cl_cruise,cd_interp_1_1)
scatter(cl_cruise,cd_interp_1_2)
grid on
xlabel("$C_L$","FontSize",14)
ylabel("$C_D$","FontSize",14)
title("Tempest/Twister Drag Polar Using Synthetic Data, EQN 1.1, and EQN 1.2 With Operating Points","FontSize",14)
legend("Synthetic","EQN 1.1","EQN 1.2","Synthetic Operating Point",...
    "EQN 1.1 Operating Point","EQN 1.2 Operating Point")


% Begin code block for 4b
% First get a power required vector
v_vec = linspace(12,60);
q_inf = 0.5.*rho_inf.*v_vec.^2;
cl_pr = W./(0.5.*rho_inf.*S.*v_vec.^2);
CD_pr0 = interp1(CL_2d,CD_2d,0);
CD_pr = CD0_CFD + (cl_pr.^2)./(pi*e*16.5);
[val23,index23] = min(CD_2d);
% CD_pr = val23 + K.*((cl_pr) -cl_pr(index23)).^2;
P_A = 0.5.*100.*26.1.*ones(1,length(v_vec));
P_R = sqrt((2.*(W.^3).*CD_pr.^2)./(rho_inf.*S.*cl_pr.^3));
% P_R = q_inf.*S.*0.025099.*v_vec + q_inf.*S.*v_vec.*(cl_pr.^2)./(pi.*e.*16.5);

figure
hold on
plot(v_vec,P_R)
plot(v_vec,P_A)
xlabel("Velocity [m/s]")
ylabel("Power [W]")
title("Tempest/Twister Power Required vs Velocity")
legend("Power Required","Power Available")
grid on

% Find and print the max velocity
[~,in3] = min(abs(P_A-P_R));
fprintf("Maximum velocity is: %f [m/s]\n\n",v_vec(in3));

% Find the max excessive power
max_excess = abs(P_A(1)-P_R(1));
RoC = max_excess/W;
fprintf("The max rate of climb is: %f [m/s]\n\n",RoC);

% Get the R over C vector and plot it against velocity
RoC = abs(P_R - P_A);

figure
hold on
plot(v_vec(1:in3),RoC(1:in3))
xlabel("Velocity [m/s]")
ylabel("Rate of Climb [m/s]")
title("Tempest/Twister Rate of Climb vs Velocity")
grid on

Rt = 1;
n = [1,1.3];
C = 18;
k = 1/(pi*e*16.5);
for i = 1:length(n)
    E(i,:) = (Rt.^(1-n(i))).*((0.5.*26.1.*C)./(q_inf.*v_vec.*S.*0.025099+...
        (2.*(W^2).*k)./(rho_inf.*v_vec.*S))).^n(i);
end

figure
hold on
grid on
plot(v_vec,E(1,:))
plot(v_vec,E(2,:))
xlabel("Velocity [m/s]")
ylabel("Endurance [h]")
title("Tempest/Twister Endurance vs Cruise Speed Range")
legend("n = 1","n = 1.3")
hold off

