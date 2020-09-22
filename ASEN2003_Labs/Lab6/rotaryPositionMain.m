%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003: Lab 6 Rotary Position Control
%   Group 13:
%   Sam Hartman
%   Joshua Seedorf
%   Johnathan Tucker
%   Sean Yoo
%   
% The purpose of this script is to model and plot the behavior of both the
% rigid and the flexible arm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear
close all
clc
%% USER INPUTS
% Gains
Kp = 50; %Proportional Gain
Kd = 1.5; %Derivative Gain
[zeta,w_n] = getZeta(Kp,Kd);
t_s = -log(0.05)./(zeta*w_n);
t_s = 0.1438;
% Target Theta
thetad = pi/4; %[rad]

%% Initial values
% Motor Components
Kg = 33.3; %Gear Ratio
Km = 0.0401; %Motor Constant [N*m/amp]
Rm = 19.2; %Armature Resistance [ohms]

% Rigid Arm Components
J_hub = 0.0005; %Base Inertia [kg*m^2]
J_load = 0.0015; %Load Inertia [kg*m^2]

% Flexible Arm Components
L = 0.45; %Length of arm [m]
M_arm = 0.06; %Mass of arm [kg]
M_tip = 0.05; %Mass on tip [kg]
fc = 1.8; %Natural frequency [Hz]

%% Initial calculations
% Inertias
J = J_hub + J_load; %Total Rigid
J_arm = (M_arm*L^2)/3; %Flexible arm
J_M = M_tip*L^2; %Tip mass
J_L = J_arm+J_M; %Total flexible link

% Stiffness
Karm = (2*pi*fc)^2*(J_L); %Flexible link stiffness

% Transfer Function Values
n1 = Kp*( (Kg*Km)/(J*Rm) ); %numerator coefficient1
d2 = 1;
d1 = ( (Kg^2*Km^2)/(J*Rm) ) + Kd*( (Kg*Km)/(J*Rm) );
d0 = Kp*( (Kg*Km)/(J*Rm) );

%% Main Simulation
% Closed loop system
num = n1;
den = [d2 d1 d0];
sysTF = tf(num,den);

% Step Response
[x,t] = step(sysTF);
x = thetad*x;
figure(4)
clf
hold on
plot(t,x)
plot([min(t) max(t)],[thetad thetad],'--')
plot([t_s t_s],[0 1.05*thetad])
plot([0 max(t)],1.05*thetad*[1 1],'k--')
plot([0 max(t)],0.95*thetad*[1 1],'k--')
sub = sprintf('\\theta_d = %.2f rad',thetad);
sub = strcat("$",sub,"$");
% title({'Step Response:',sub});
% ylabel('\theta [rad]');
ylabel('$\theta\:[rad]$', 'interpreter','latex','FontSize',20)
xlabel('$Time\:[s]$', 'interpreter','latex','FontSize',20)
title({'$Arm\:Angular\:Position\:Over\:Time\:(Modeled)$',sub},...
    'Interpreter','latex','FontSize',20)
legend("Modeled Behavior","Desired Theta","Settling Time","5% Margin")

%% Flex plots
thetad = 0.1;
% Flex constants
Kg = 33.3; %Gear Ratio
Km = 0.0401; %Motor Constant [N*m/amp]
Jhub = 0.0005; %Base Inertia [kg*m^2]
L = 0.45; %Length of arm [m]
Mtip=0.05;
JM=Mtip*L^2;
Marm=0.06;
Jarm=(Marm*L^2)/3;
JL = Jarm+JM; %L
Rm = 19.2; %Armature Resistance [ohms]
fc = 1.8; %Natural frequency [Hz]
Karm = (2*pi*fc)^2*(JL); %Flexible link stiffness

% Create constants for transfer function denominator
p1= -(Kg^2 * Km^2) / (Jhub *Rm);
p2= (Kg^2 * Km^2 *L) / (Jhub *Rm);
q1=Karm / (L*Jhub);
q2= - (Karm*(Jhub+JL)) / (JL*Jhub);
r1=(Kg*Km) / (Jhub * Rm);
r2= -( Kg*Km*L) / (Jhub *Rm);

%Random values
K1= 15.5;   %Kptheta;
K2= 0;  %Kpdisplace; 
K3= 1.5;  %KDtheta;
K4= 1; %KDdisplace;

[zeta,w_n] = getZeta(15.5,1.5);
t_s = -log(0.05)./(zeta*w_n);
t_s = 0.3767;
% t_s = 0.2830;
% t_s = 0.2450;

% Create coefficients for transfer function denominator
lambda3=-p1+K3*r1+K4*r2;
lambda1=p1*q2-q1*p2+K3*(q1*r2-r1*q2)+K2*(p2*r1-r2*p1);
lambda2=-q2+K1*r1+K2*r2+K4*(p2*r1-r2*p1);
lambda0=K1*(q1*r2-r1*q2);

% Create coefficients for transfer function numerator
cn1_2=K1 * r1;
cn1_0=K1*(q1*r2-r1*q2);
cn2_2=K1 * r2;
cn2_1=K1*(p2*r1 -r2*p1);

% Create constants for transfer function denominator
cd4=1;
cd3=lambda3;
cd2=lambda2;
cd1=lambda1;
cd0=lambda0;

%% Main Simulation
% Closed loop system
num = [cn1_2 0 cn1_0];
den = [cd4 cd3 cd2 cd1 cd0];
sysTF = tf(num,den);
% Step Response
[x,t] = step(sysTF);
x = thetad*x;
figure(1)
clf
hold on
plot(t,x)
plot([min(t) max(t)],[thetad thetad],'--')
plot([t_s t_s],[0 1.1*thetad])
plot([0 max(t)],1.1*thetad*[1 1],'k--')
plot([0 max(t)],0.90*thetad*[1 1],'k--')
sub = sprintf('\\theta_d = %.2f rad',thetad);
sub = strcat("$",sub,"$");
% title({'Step Response:',sub});
% ylabel('\theta [rad]');
ylabel('$\theta\:[rad]$', 'interpreter','latex','FontSize',20)
xlabel('$Time\:[s]$', 'interpreter','latex','FontSize',20)
title({'$Arm\:Angular\:Position\:Over\:Time\:(Modeled)$',sub},...
    'Interpreter','latex','FontSize',20)
legend("Modeled Behavior","Desired Theta","Settling Time","10% Margin")

%% Tip deflection plot
num = [0 cn2_2 cn2_1];
den = [cd4 cd3 cd2 cd1 cd0];
sysTF = tf(num,den);
% Step Response
[x,t] = step(sysTF);
x = thetad*x;

figure(3)
plot(t,x)
% plot([min(t) max(t)],[thetad thetad],'--')
% plot([t_s t_s],[0 1.1*thetad])
% plot([0 max(t)],1.1*thetad*[1 1],'k--')
% plot([0 max(t)],0.90*thetad*[1 1],'k--')
sub = sprintf('\\theta_d = %.2f rad',thetad);
sub = strcat("$",sub,"$");
% title({'Step Response:',sub});
% ylabel('\theta [rad]');
xlim([0 1]);
ylabel('$Tip\:Deflection\:[m]$', 'interpreter','latex','FontSize',20)
xlabel('$Time\:[s]$', 'interpreter','latex','FontSize',20)
title({'$Arm\:Tip\:Deflection\:Over\:Time\:(Modeled)$',sub},...
    'Interpreter','latex','FontSize',20)
legend("Modeled Behavior","Desired Theta","Settling Time","10% Margin")