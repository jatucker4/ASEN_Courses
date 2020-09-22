%% 


% Housekeeping
clear all;

% Variable declaration

mu_Su = 132712000000; %(km^3/s^2)

mu_E = 398600; %(km^3/s^2)
a_E = 149.6 *10^6; %km

mu_Ma = 37931000; %(km^3/s^2)
a_Ma = 227.9 *10^6; %km
a_t = .5 * (a_E + a_Ma); %km


h_E = sqrt(mu_Su * a_E); %km^2/s
h_Ma = sqrt(mu_Su * a_Ma); %km^2/s
h_t = sqrt(2 * mu_Su) * sqrt((a_E*a_Ma)/(a_E + a_Ma));

v_E = h_E/a_E; %km/s
v_Ma = h_Ma/a_Ma; %km/s
v_t_E = h_t/a_E;
v_t_Ma = h_t/a_Ma;

delv_t_1 = abs(v_t_E - v_E);
delv_t_2 = abs(v_t_Ma - v_Ma);
delv_t_tot = delv_t_1 + delv_t_2;

t_t = .5*((2*pi)/sqrt(mu_Su)) * a_t^(3/2);
T_Earth = 2*pi*sqrt((a_E^3) / mu_Su);
T_Mars = 2*pi*sqrt((a_Ma^3) / mu_Su);
t_transfer_total = t_t + T_Earth + T_Mars;




