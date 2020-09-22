%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003: Lab 6 Rotary Position Control
%   Group 13:
%   Sam Hartman
%   Joshua Seedorf
%   Johnathan Tucker
%   Sean Yoo
%   
% The purpose of this script is to run a Monte Carlo simulation to 
% determine the optimal flexible gains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
close all; clear all; clc;
%% Create constants
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

p1= -(Kg^2 * Km^2) / (Jhub *Rm);
p2= (Kg^2 * Km^2 *L) / (Jhub *Rm);
q1=Karm / (L*Jhub);
q2= - (Karm*(Jhub+JL)) / (JL*Jhub);
r1=(Kg*Km) / (Jhub * Rm);
r2= -( Kg*Km*L) / (Jhub *Rm);

%Random values for the simulation
K1_vec= linspace(0,20,20);   %Kptheta;
K2_vec= linspace(-50,0,20);   %Kpdisplace; 
K3_vec= linspace(0,1.5,20);  %KDtheta;
K4_vec= linspace(0,1.5,20);  %KDdisplace;

thetaD=0.1;
%Loop through every vector value
for i = 1:10
    K1 = K1_vec(i);
    for j = 1:10
        K2 = K2_vec(j);
        for k = 1:10
            K3 = K3_vec(k);
            K4 = K4_vec(k);
            % Create denominator constants
            lambda3=-p1+K3*r1+K4*r2;
            lambda1=p1*q2-q1*p2+K3*(q1*r2-r1*q2)+K2*(p2*r1-r2*p1);
            lambda2=-q2+K1*r1+K2*r2+K4*(p2*r1-r2*p1);
            lambda0=K1*(q1*r2-r1*q2);
            % Create numerator constants
            cn1_2=K1 * r1;
            cn1_0=K1*(q1*r2-r1*q2);

            cn2_2=K1 * r2;
            cn2_1=K1*(p2*r1 -r2*p1);

            
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
            x = thetaD*x;
            % Find if any values surpass the time threshold
            time_threshold = find(abs(x-thetaD) < thetaD * .1);
            % If no values surpass time threshold then add it to viable
            % pairs matrix
            if ~isempty(time_threshold)
                pair{j,k} = [K1,K2,K3,K4];
            end
        end
    end
end