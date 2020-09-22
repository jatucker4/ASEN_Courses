function [zeta,w_n] = getZeta(Kp, Kd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Created by: Johnathan Tucker
%
% The purpose of this function is to calculate and output the zeta and
% omega n values given the proportional and derivative gains
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define constants
Kg = 33.3; %Gear Ratio
Km = 0.0401; %Motor Constant [N*m/amp]
Rm = 19.2; %Armature Resistance [ohms]
% Rigid Arm Components
J_hub = 0.0005; %Base Inertia [kg*m^2]
J_load = 0.0015; %Load Inertia [kg*m^2]
J = J_hub + J_load; %Total Rigid
%% Calculate Zeta
zeta = (Kg.^2 .* Km^2  + (Kd .* Kg .* Km))./...
    (2.*sqrt(Kp .* Kg .* Km .* J .* Rm));
%% Omega N
w_n = sqrt((Kg .* Km .* Kp)./(J .*Rm));
end

