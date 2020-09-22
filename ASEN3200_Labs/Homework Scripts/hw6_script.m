%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is for the 6th homework assignment in ASEN 3200
%
% Created by: Johnathan Tucker
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
close all;
clear all;
%% Create vectors and matrix for Q5 part a
% Get period
T = 2*pi*sqrt((6878^3)/398600.4415);
% Get mean motion
n = 2*pi/T;
% Create CW matrices
phi_rr_test = [(4 - 3*cos(pi/2)) , 0 , 0;6*(sin(pi/2) - pi/2) , 1, 0;...
        0, 0, cos(pi/2)];
phi_rv_test = [(1/n)*sin(pi/2) , (2/n)*(1-cos(pi/2)), 0;...
        (2/n)*(cos(pi/2) - 1), (1/n)*(4*sin(pi/2) - 3*pi/2), 0;...
        0, 0, (1/n)*sin(pi/2)];
% Create r_0
temp1 = [0.25; -0.4; 0];
% Create first multiplication solution
temp2 = phi_rr_test*temp1;
% Create desired position vector
b = [0; -0.1; 0] - temp2;
% Solve for velocity
d_v_plus_0 = phi_rv_test\b;
% Check to make sure it works 
test = phi_rr_test * temp1 + phi_rv_test*d_v_plus_0;
d_v_minus_0 = [0 ; -0.1e-3; 0 ];
delta_v_1 = d_v_plus_0 - d_v_minus_0;

%% Get the velocity at T/4 for part b and then delta v 2
d_v_minus_f = [3*n*.250; 6*n*.4; 0] + [-2*d_v_plus_0(2); (-2*d_v_plus_0(1)...
    - 3*d_v_plus_0(2)); 0];

%% Create the code for the plots in part c
% Linearly spaced time vector
t = linspace(0,T/4);
r_0 = [0.250 ; -0.4 ; 0];
% Loop through each time values solving CW eqns analytically each time
for i = 1:length(t)
    % Create CW matrices
    phi_rr = [(4 - 3*cos(n*t(i))) , 0 , 0;6*(sin(n*t(i)) - n*t(i)) , 1, 0;...
        0, 0, cos(n*t(i))];
    phi_rv = [(1/n)*sin(n*t(i)) , (2/n)*(1-cos(n*t(i))), 0;...
        (2/n)*(cos(n*t(i)) - 1), (1/n)*(4*sin(n*t(i)) - 3*n*t(i)), 0;...
        0, 0, (1/n)*sin(n*t(i))];
    % Get the position at time t
    r_t(:,i) = phi_rr*r_0 + phi_rv*d_v_plus_0;
    
    % Create the rest of the CW matrices
    phi_vr = [3*n*sin(n*t(i)), 0, 0; (6*n)*(cos(n*t(i)) -1), 0, 0;...
        0, 0, -n*sin(n*t(i))];
    phi_vv = [cos(n*t(i)), 2*sin(n*t(i)), 0;...
        -2*sin(n*t(i)), 4*cos(n*t(i))-3, 0; 0, 0, cos(n*t(i))];
    v_t(:,i) = phi_vr*r_0 + phi_vv*d_v_plus_0;
end

figure(1)
plot(r_t(1,:),r_t(2,:))
hold on
plot(0,0,'o')
xlabel("$X\:position[km]$",'Interpreter','latex','FontSize',26)
ylabel("$Y\:position[km]$",'Interpreter','latex','FontSize',26)
legend("Relative B Position", "Relative A Position")
title("$XY-Position\:of\:CubeSat\:B\:Relative\:to\:A\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)

figure(2)
plot(v_t(1,:),v_t(2,:))
hold on
plot(0,0,'o')
xlabel("$X\:velocity[km/s]$",'Interpreter','latex','FontSize',26)
ylabel("$Y\:velocity[km/s]$",'Interpreter','latex','FontSize',26)
legend("Relative B Velocity", "Relative A Velocity")
title("$XY-Velocity\:of\:CubeSat\:B\:Relative\:to\:A\:Over\:Time$",...
    'Interpreter','latex','FontSize',26)