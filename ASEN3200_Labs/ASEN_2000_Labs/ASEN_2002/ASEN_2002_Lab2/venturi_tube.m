function velocity = venturi_tube(pressure_change, T_atm, P_atm)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

pressure_change = abs(pressure_change);
velocity = sqrt((2*pressure_change*287*T_atm)/(P_atm*(1-(1/9.5)^2)));
end

