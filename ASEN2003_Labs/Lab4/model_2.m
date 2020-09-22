function omega = model_2(moment,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses work-energy to derive the model equation for the 
% unbalanced wheel apparatus also including the negative moment due to
% friction.
%
% Created by: Johnathan Tucker
%
% Inputs:
%           moment: The negative moment due to friction (N-m)
%           theta: Angluar position of the wheel (radians)
%
% Outputs:
%           omega: The angular velocity of the wheel with the particle
%                   (radians/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Constants
M = 11.7; % mass of cylinder (kg)
M0 = 0.7; % mass of trailing supports (kg)
m = 3.4; % mass of extra mass (kg)
R = 0.235; % radius of cylinder (m)
k = 0.203; % radius of gyration of wheel (m)
I = M*k^2; % moment of Inertia
beta = 5.14*(pi/180); % slope of ramp (rad)
r = 0.178; % radius to extra mass (m)
rem = 0.019; % radius of extra mass (m)
g = 9.81; % gravitational acceleration (m/s^2)
h1 = 0.37338 + R; % start height of of cyclinder (m)
%% Calculate the omega value using the formula from the lab derivation
omega = sqrt((2.*theta.*(g.*R.*sin(beta).*(M+M0)-moment))/(I + R^2 * (M+M0)));
end