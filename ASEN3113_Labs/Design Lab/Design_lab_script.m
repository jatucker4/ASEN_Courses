%% Housekeeping
clc;
clear all;
close all;
% Create a model for the earth's distance from the sun
%% Define constants
boltz = 5.67e-8;
a = 149.6e6; % [km]
mu = 132712000000;% [km]
e = 0.0167;
x_0 = a*(1-e); % [km]
y_0 = 0; % [km]
z_0 = 0; % [km]
xdot_0 = 0; % [km/s]
ydot_0 = sqrt(2*((mu/x_0) - (mu/(2*a)))); % [km/s]
zdot_0 = 0; % [km/s]
T = 2*pi*sqrt((a^3)/mu);
state_init = [x_0, y_0, z_0, xdot_0, ydot_0, zdot_0];
tspan = [0 T];
%% Plug into ode45
[t,y] = ode45(@(t,y) odefunc(t,y), tspan, state_init);
%% Verify via plot
figure(1)
plot(y(:,1),y(:,2))
hold on
scatter(0,0)
title("Earth's Orbit About the Sun")
legend("Earth Orbit", "Sun")
xlabel("x-axis (km)")
ylabel("y-axis (km)")
zlabel("z-axis (km)")
%% Get a vector of positions
r = sqrt(y(:,1).^2 + y(:,2).^2);
% Convert to Au
r = r./1.496e8;
%% Get radiation
G = 1361./(r.^2);
%% Get Ir backload
Wrate = 88;
Srate = 63;
Wd = 1-e;
Sd = 1+e;

%% Get the max q_in
days = t./86400;
% First get a vector of angles that the precession goes through
angle = 23.45*sin((360/365)*(284 + days));
a = 0.2;
epsilon = 0.85;
for i = 1:length(r)
    ir_backload = ((Wrate - Srate)/(Wd-Sd))*(r(i)-Sd) + Srate;
    ir_backload = 75.5;
    G_temp = 1407.622;
%     Q_in(i) = ir_backload*epsilon + G(i)*a*cosd(angle(i));
end
ang_perc = [linspace(0,100),linspace(100,0)]/100;
Q_in = ir_backload*epsilon + G_temp*a*cosd(23.5)*(ang_perc);
[Max, distance] = max(Q_in);

fprintf("Max radiation is: %f W/m^2\n",Max);
fprintf("Max radiation occurs at: %f Au\n",r(distance));

%% Create a symbolic function for A
% syms sig G_s
% A = 20/(boltz*sig*T^4 - ir_backload*sig - a*G_s*cosd(angle));
T = 30 + 273.15;
A = 20/(boltz*epsilon*T^4 - ir_backload*epsilon - a*G_temp*cosd(23.5));


time = 0:1:12*3600;
time_equinox = [0:1:11*3600 , zeros(1,3600)];

theta = 360/86400;

Gs_winter = 1361/(1-0.0167)^2;
Gs_Summer = 1361/(1+0.0167)^2;
Gs_Equinox = 1361;

ir_equinox = mean([88,63]);
A = 20/(boltz*epsilon*T^4 - 75.5*epsilon - a*Gs_Equinox*cosd(0));

Gs1 = sind(theta*time)*Gs_winter;
Gs2 = sind(theta*time)*Gs_Summer;
Gs3 = sind(theta*time_equinox)*Gs_Equinox;
time = [time,length(time):1:2*length(time)-1];
Gs1 = [Gs1,zeros(1,length(Gs1))];
Gs2 = [Gs2,zeros(1,length(Gs2))];
Gs3 = [Gs3,zeros(1,length(Gs3))];

figure(2)
plot(time,Gs1)
hold on
plot(time,Gs2)
hold on
plot(time,Gs3)
xticks([0 length(time_equinox) length(Gs1)])
xticklabels({'Noon','Midnight','Noon'})



temp_vec = ((20 + A*75.5*epsilon + A.*Gs3.*a)./(A*epsilon*boltz)).^(1/4);% double(solve(A*epsilon*(T^4) == 20 + A*75.5 + A*Gs3*a));
power_drawn = A.*epsilon.*boltz.*(temp_vec.^4);
power_drawn = 20 + A*75.5*epsilon + A.*Gs3.*a;
temp_vec = temp_vec - 273.15;

figure(3)
plot(time,temp_vec);
hold on
plot(time,power_drawn)
xticks([0 length(time_equinox) length(Gs1)])
xticklabels({'Noon','Midnight','Noon'})
