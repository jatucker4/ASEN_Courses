%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003: Lab 6 Rotary Position Control
%   Group 13:
%   Sam Hartman
%   Joshua Seedorf
%   Johnathan Tucker
%   Sean Yoo
%   
% The purpose of this script is to run a Monte Carlo simulation to
% determine the viable gains for the rigid arm system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
close all;
clear all;
%% Define constants
Kg = 33.3; %Gear Ratio
Km = 0.0401; %Motor Constant [N*m/amp]
Rm = 19.2; %Armature Resistance [ohms]
% Rigid Arm Components
J_hub = 0.0005; %Base Inertia [kg*m^2]
J_load = 0.0015; %Load Inertia [kg*m^2]
J = J_hub + J_load; %Total Rigid
volt_constraint = 10;
time_constraint = 0.15;
N = 100;
experiment = "rigid"; 
thetaD = pi/4;
%% Start simulation
if contains(experiment,"rigid")
    Kp_max = 50;
    Kp_min = -2;
    Kp_vec = linspace(Kp_min,Kp_max,N);
    
    Kd_max = 1.5;
    Kd_min = -1.5;
    Kd_vec = linspace(Kd_min,Kd_max,N);
    
    for i = 1:N
        Kp = Kp_vec(i);
        for j = 1:N
            Kd = Kd_vec(j);
            % Transfer Function Values
            n1 = Kp*( (Kg*Km)/(J*Rm) ); %numerator coefficient1
            d2 = 1;
            d1 = ( (Kg^2*Km^2)/(J*Rm) ) + Kd*( (Kg*Km)/(J*Rm) );
            d0 = Kp*( (Kg*Km)/(J*Rm) );
            % Closed loop system
            num = n1;
            den = [d2 d1 d0];
            sysTF = tf(num,den);
            % Step Response
            [x,t] = step(sysTF);
            x = thetaD*x;
            % Get the v_in vector
            temp_1 = Kp .* (thetaD - x);
            temp_1(1) = [];
            temp_2 = (Kp .* (diff(x)./diff(t)));
            V_in = temp_1 - temp_2;
            % Check if any values surpass v_in
            time_threshold = find(abs(x-thetaD) < thetaD * .05);
            v_limit_check = find(abs(V_in) > volt_constraint);
            v_deadband_check = find(abs(V_in) > 1.5);
            bool_test = ~isempty(time_threshold);
            % Test if the Kp and Kd values surpass limits
            if bool_test && isempty(v_limit_check) && isempty(v_deadband_check)
                % If not add them to the list of viable gains
                    reach_rigid = t(time_threshold(1));
                    viable_pairs{i,j} = [Kd,Kp];
                    viable_Kd(i,j) = Kd;
                    viable_Kp(i,j) = Kp;
            end
        end
    end
end  



