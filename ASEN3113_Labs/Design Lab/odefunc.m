function dydt = odefunc(t,y)
%This is the odefunction for the first lab of ASEN3200
% Create output vector
dydt = zeros(6,1);
% Create gravitational parameter
mu = 132712000000;
% Get the position for the orbit
r = sqrt(y(1)^2 + y(2)^2 + y(3)^2);
% Create the outputs
% Note: the first three outputs are just the velocities and the last three
% are taken from page 69 in the book and are the defined accelerations in
% each direction
dydt(1) = y(4);
dydt(2) = y(5);
dydt(3) = y(6);
dydt(4) = -(mu/r^3) * y(1);
dydt(5) = -(mu/r^3) * y(2);
dydt(6) = -(mu/r^3) * y(3);
end
