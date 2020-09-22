function dydt = odefunc(t,y,wind,m)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
dydt = zeros(6,1);
rho = 1.225; % density of air at sea level standard atmosphere
S = (pi*.03^2)/4;
Cd = 0.6;
g = [0  0 -9.81];
% First get the relative wind by adding the inertial and wind velocity
% vectors
v_inertial = [y(4), y(5), y(6)];
v = v_inertial - wind;

% Now calculate the drag force
F_d = 0.5*rho*S*Cd*norm(v)^2;
F_d = -(v/norm(v))*F_d;

% Now calculate force due to gravity
F_g = g*m;

% Sum the forces
F = F_g + F_d;

% Get the accelerations
a_vec = F/m;

% Output the state vector
dydt(1) = y(4);
dydt(2) = y(5);
dydt(3) = y(6);
dydt(4) = a_vec(1);
dydt(5) = a_vec(2);
dydt(6) = a_vec(3);
end

