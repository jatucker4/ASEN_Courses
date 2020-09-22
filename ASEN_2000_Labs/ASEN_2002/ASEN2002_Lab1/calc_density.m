function density = calc_density(radius,F_System,mass_Helium)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Define constants
sigma_SB = 5.670*10^(-8);
alpha_SB = 0.6;
epsilon_B = 0.8;
alpha_EB = epsilon_B;
q_sun = 1353;
q_earth = 237;

%% Solving for Temperature with Radiation
Q_solar = alpha_SB*q_sun*pi*radius^2;
Q_earth = alpha_EB*q_earth*pi*radius^2;

Temp = ((Q_solar + Q_earth)/(epsilon_B*sigma_SB*4*pi*radius^2))^(1/4);

%% Solving for Volume
Volume = (mass_Helium*2.0769*Temp)/0.5589235;

%% Solve for density
density = F_System/(Volume*9.81);
% From this density we have to do a reverse look up in atmoscoesa where
% We guess and check to find the altitude that matches this density
end

