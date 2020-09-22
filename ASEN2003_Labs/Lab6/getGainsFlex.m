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
            zeta = getZeta(Kp,Kd);
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

            % Check if voltage is within the limit
            temp_1 = Kp .* (thetaD - x);
            temp_1(1) = [];
            temp_2 = (Kp .* (diff(x)./diff(t)));
%             V_in = Kp .* (thetaD - x) - (Kp .* (diff(x)./diff(t)));
            V_in = temp_1 - temp_2;
            time_threshold = find(abs(x-thetaD) < thetaD * .05);
            v_limit_check = find(abs(V_in) > volt_constraint);
            bool_test = ~isempty(time_threshold);
            if bool_test
%                 for k = 1:length(V_in)
                    reach_rigid = t(time_threshold(1));
                    viable_pairs{i,j} = [Kd,Kp];
                    viable_Kd(i,j) = Kd;
                    viable_Kp(i,j) = Kp;
%                 end
%             elseif max(abs(V_in)) > 10.05 && max(abs(V_in)) < 10.05
%                 for k = 1:length(V_in)
%                     if bool_test
%                         viable_Kd(i,j) = Kd;
%                         viable_Kp(i,j) = Kp;
%                     end
%                 end
            end
        end
    end
end  