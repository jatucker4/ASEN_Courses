%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 9 - Integrating Forward in Time with ODE45
% 
% The purpose of this challenge is to solve for a car's velocity with
% respect to time. We will use the built-n ODE 45 function to solve
% over 30 seconds with varying accelerations. Three parts will include
% accelerating with gas, N2O, and braking until the car comes to a stop.
%
% This code should also incorporate your previous integration methods from
% the previous coding challenges, so you can compare side-by-side the
% accuracies of the methods.
%
% Please ZIP and upload your team's script(s) and figures to Canvas to 
% complete the challenge.
% 
% STUDENT TEAMMATES
% 1. Johnathan Tucker
% 2. Daniel Denton
% 3.
% 4.
%
% CHALLENGE AUTHORS
% Torin Clark, Melinda Zavala, Justin Fay, John Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clearvars
close all
clc

%% Given Variables

Cd   = 0.5;   % coefficient of drag
rho  = 1.225; % density of air at sea level, at 15C (kg/m^3)
aGas = 3;     % acceleration with gas (m/s^2)
aN2O = 5;     % acceleration with NO2 (m/s^2)
aB   = -1;    % braking acceleration (m/s^2)
Ar   = 4;     % frontal area of vehicle (m^2)
m    = 800;   % kg

%% Part 1: Acceleration
% Solve for the acceleration with gas for the first 15 seconds
% Function: dv/dt = a(t) - (1/2)*(1/m)*Cd*rho*A*v(t)^2 (second component is
% equation for drag)

% IMPORTANT:
% Your code should run 4th order RK along with Euler's method, so you can
% compare the two. Use your Euler method code from the previous coding
% challenge to avoid re-writting Euler's method.
h = 0.1;
v(1,1) = 0;
time = 0:.1:30;
a_part1 = 3;
for i = 2:length(time)
    v(1,i) = v(1,i-1) + h*((a_part1 - (1/2)*(Ar/m)*Cd*rho*v(1,i-1)^2));
end
% Initialize acceleration function for only gas
a_gas = @(t,v) aGas - (1/2)*(1./m).*Cd.*rho.*Ar.*v.^2;
v_r_1(1,1) = 0;
time_r_1 = 0:0.1:30;
for i = 2:length(time_r_1)
    k_1 = a_gas(time_r_1(1,i-1),v_r_1(1,i-1));
    k_2 = a_gas(time_r_1(1,i-1) + h/2, v_r_1(1,i-1) + k_1*h/2);
    k_3 = a_gas(time_r_1(1,i-1) + h/2, v_r_1(1,i-1) + k_2*h/2);
    k_4 = a_gas(time_r_1(1,i-1) + h, v_r_1(1,i-1) + k_3*h);
    v_r_1(1,i) = v_r_1(1,i-1) + (h/6)*(k_1 +2*k_2 + 2*k_3 +k_4);
end
tspan_1 = [0 30];
[t_1, y_1] = ode45(a_gas, tspan_1, 0);
figure(1)
plot(time, v)
hold on
plot(time_r_1, v_r_1)
hold on
plot(t_1,y_1)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity Versus Time')
legend('Euler', 'Runge Kutta', 'ode45')
%% Part 2: N20
% If press N20 pedal from t = 15-18 seconds. Plot until t = 18.
v2(1,1) = 0;
a_part2 = 5;
time_part2 = 0:.1:18;
for i = 2:length(time)
    if time(i) < 15
        v2(1,i) = v2(1,i-1) + h*((a_part1 - (1/2)*(Ar/m)*Cd*rho*v2(1,i-1)^2));
    else
        v2(1,i) = v2(1,i-1) + h*((a_part2 - (1/2)*(Ar/m)*Cd*rho*v2(1,i-1)^2));
    end
end
% Initialize acceleration function for N2O
a_N2O = @(t,v) aN2O - (1/2)*(1./m).*Cd.*rho.*Ar.*v.^2;
v_r_2(1,1) = 0;
time_r_2 = 0:0.1:30;
for i = 2:length(time_r_2)
    if time_r_2(1,i) < 15
        k_1 = a_gas(time_r_2(1,i-1),v_r_2(1,i-1));
        k_2 = a_gas(time_r_2(1,i-1) + h/2, v_r_2(1,i-1) + k_1*h/2);
        k_3 = a_gas(time_r_2(1,i-1) + h/2, v_r_2(1,i-1) + k_2*h/2);
        k_4 = a_gas(time_r_2(1,i-1) + h, v_r_2(1,i-1) + k_3*h);
        v_r_2(1,i) = v_r_2(1,i-1) + (h/6)*(k_1 +2*k_2 + 2*k_3 +k_4);
    else
        k_1 = a_N2O(time_r_2(1,i-1),v_r_2(1,i-1));
        k_2 = a_N2O(time_r_2(1,i-1) + h/2, v_r_2(1,i-1) + k_1*h/2);
        k_3 = a_N2O(time_r_2(1,i-1) + h/2, v_r_2(1,i-1) + k_2*h/2);
        k_4 = a_N2O(time_r_2(1,i-1) + h, v_r_2(1,i-1) + k_3*h);
        v_r_2(1,i) = v_r_2(1,i-1) + (h/6)*(k_1 +2*k_2 + 2*k_3 +k_4);
    end
end
tspan2 = [0 18];
[t_2,y_2] = ode45(a_N2O, tspan2, 0);
figure(2)
plot(time_r_2(1,1:183),v_r_2(1,1:183))
hold on
plot(time_part2(1,1:181),v2(1,1:181))
hold on
plot(t_2,y_2)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity Versus Time with NO2')
legend('Runge Kutta', 'Eulers','ode45')

%% Part 3: deceleration
% Use accelertion value of -1 m/s^2, plot until v = 0 (could be > 30sec)
a_part3 = -1;
time_2 = 0:.1:100;
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
% Initialize acceleration function for braking
a_brake = @(t,v) (aB - (1/2)*(1./m).*Cd.*rho.*Ar.*v.^2).*(v>0);
time_r_2 = 0:0.1:100;
for i = 2:length(time_r_2)
    if time_r_2(1,i) < 15
        k_1 = a_gas(time_r_2(1,i-1),v_r_2(1,i-1));
        k_2 = a_gas(time_r_2(1,i-1) + h/2, v_r_2(1,i-1) + k_1*h/2);
        k_3 = a_gas(time_r_2(1,i-1) + h/2, v_r_2(1,i-1) + k_2*h/2);
        k_4 = a_gas(time_r_2(1,i-1) + h, v_r_2(1,i-1) + k_3*h);
        v_r_2(1,i) = v_r_2(1,i-1) + (h/6)*(k_1 +2*k_2 + 2*k_3 +k_4);
    elseif time_r_2(1,i) > 14 && time_r_2(1,i) < 19
        k_1 = a_N2O(time_r_2(1,i-1),v_r_2(1,i-1));
        k_2 = a_N2O(time_r_2(1,i-1) + h/2, v_r_2(1,i-1) + k_1*h/2);
        k_3 = a_N2O(time_r_2(1,i-1) + h/2, v_r_2(1,i-1) + k_2*h/2);
        k_4 = a_N2O(time_r_2(1,i-1) + h, v_r_2(1,i-1) + k_3*h);
        v_r_2(1,i) = v_r_2(1,i-1) + (h/6)*(k_1 +2*k_2 + 2*k_3 +k_4);
    else
        k_1 = a_brake(time_r_2(1,i-1),v_r_2(1,i-1));
        k_2 = a_brake(time_r_2(1,i-1) + h/2, v_r_2(1,i-1) + k_1*h/2);
        k_3 = a_brake(time_r_2(1,i-1) + h/2, v_r_2(1,i-1) + k_2*h/2);
        k_4 = a_brake(time_r_2(1,i-1) + h, v_r_2(1,i-1) + k_3*h);
        v_r_2(1,i) = v_r_2(1,i-1) + (h/6)*(k_1 +2*k_2 + 2*k_3 +k_4);
    end
end
tspan3 = [19 100];
[t_1, y_1] = ode45(a_gas, [0 14], 0);
[t_2,y_2] = ode45(a_N2O, [15 18], y_1(end,1));
[t_3,y_3] = ode45(a_brake, tspan3, y_2(end,1));
t_total = [t_1 ; t_2; t_3];
y_total = [y_1; y_2; y_3];
figure(3)
plot(time_2,v3)
hold on
plot(time_r_2, v_r_2)
hold on
plot(t_total,y_total)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity Versus Time with NO2 and Deceleration')
legend('Eulers', 'Runge Kutta', 'ode45')

%% Bonus: New Geo Metro
% You are an engineering intern working on a new low-drag profile for
% the next generation Geo Metro. You took a scale model of the car, put it
% in a wind tunnel, and collected data to determine the drag acting on the
% car at different speeds.

% This Geo Metro has a dynamic profile to reduce drag at high speeds.

% (1) Compute a fit for the drag vs. velocity data. You may use whichever
% formulation (linear, quadratic, etc.) that you think is appropriate.

% (2) Create an anonymous function to calculate the acceleration of the
% vehicle at a given velocity V.

% (3) Use ODE 45 to determine how long it takes the new Geo Metro to go
% from 0 to 60 miles per hour (mph)

% (4) Using the best estimate for the braking acceleration, how long would
% it take the Geo Metro to stop from 120 mph?

% (5) How long would it take the Geo Metro to come to a complete stop using
% air resistance alone? Does the car ever stop only due to air resistance
% with your model of drag?

ac_geo = 2; % m/s^2
ac_brake = [-0.75 -0.8 -0.79 -0.5 -0.71]; % Results from different trials

load('aeromats')
figure(4)
plot(v_test, drag_noisy, '*');
xlabel('Air Velocity (s)')
ylabel('Measured Drag (m/s$^2$)')
grid on;
title('Air Drag vs. Air Velocity')
hold on;





