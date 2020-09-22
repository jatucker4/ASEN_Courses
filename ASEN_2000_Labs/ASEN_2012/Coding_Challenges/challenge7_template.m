%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 7 - Integrating Forward in Time with Euler (RK1)
% 
% The purpose of this challenge is to solve for a car's velocity with
% respect to time. We will use Euler's Method solving forward in time (RK1)
% over 30 seconds with varying accelerations. Three parts will include
% accelerating with gas, N2O, and braking until the car comes to a stop.
%
% Please ZIP and upload your team's script(s) and figures to Canvas to 
% complete the challenge.
% 
% STUDENT TEAMMATES
% 1. Johnathan Tucker
% 2. Emerson Beinhauer
% 3. Parker Simmons
% 4. Cameron Turman
%
% CHALLENGE AUTHORS
% Torin Clark, Melinda Zavala 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all
close all
clc

%% Given Variables

Cd   = 0.5;   % coefficient of drag
rho  = 1.225; % density of air at sea level, at 15C (kg/m^3)
aGas = 3;     % acceleration with gas (m/s^2)
aN20 = 5;     % acceleration with NO2 (m/s^2)
aB   = -1;    % braking acceleration (m/s^2)
Ar   = 4;     % frontal area of vehicle (m^2)
m    = 800;   % kg
time = 0:1:30;
a_part1 = 3;
h = 0.1;
%% Part 1: Acceleration
% Solve for the acceleration with gas for the full time period and plot the
% results.
% Function: dv/dt = a(t) - (1/2)*(A/m)*Cd*rho*v(t)^2 (second component is
% equation for drag)
v(1,1) = 0;
for i = 2:length(time)
    v(1,i) = v(1,i-1) + h*((a_part1 - (1/2)*(Ar/m)*Cd*rho*v(1,i-1)^2));
end
figure(1)
plot(time,v)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity Versus Time')
%% Part 2: N20
% If press N20 pedal from t = 15-18 seconds. Plot until t = 18.
v2(1,1) = 0;
a_part2 = 5;
time_part2 = 0:1:18;
for i = 2:length(time)
    if time(i) < 15
        v2(1,i) = v2(1,i-1) + h*((a_part1 - (1/2)*(Ar/m)*Cd*rho*v2(1,i-1)^2));
    else
        v2(1,i) = v2(1,i-1) + h*((a_part2 - (1/2)*(Ar/m)*Cd*rho*v2(1,i-1)^2));
    end
end
figure(2)
plot(time_part2,v2(1,1:1:19))
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity Versus Time with NO2')
%% Part 3: deceleration
% Use accelertion value of -1 m/s^2, plot until v = 0 (could be > 30sec)
a_part3 = -1;
time_2 = 0:1:100;
v3(1,1) = 0;
for i = 2:length(time_2)
    if time_2(i) < 15
        v3(1,i) = v3(1,i-1) + h*((a_part1 - (1/2)*(Ar/m)*Cd*rho*v3(1,i-1)^2));
    elseif time_2(i) > 14 && time_2(i) < 19
        v3(1,i) = v3(1,i-1) + h*((a_part2 - (1/2)*(Ar/m)*Cd*rho*v3(1,i-1)^2));
    else
        v3(1,i) = v3(1,i-1) + h*((a_part3 - (1/2)*(Ar/m)*Cd*rho*v3(1,i-1)^2));
    end
end
end_index = find(v3<0);
v3 = v3(1, 1:end_index(1,1));
time_2 = time_2(1,1:end_index(1,1));
figure(3)
plot(time_2,v3)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity Versus Time with NO2 and Deceleration')