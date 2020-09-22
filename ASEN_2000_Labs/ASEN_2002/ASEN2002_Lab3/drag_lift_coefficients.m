function [C_d,C_l] = drag_lift_coefficients(x_pos, y_pos,C_p,aoa)
% drag_lift_coefficients
% This function takes in the x position, y position, angle of attack, and
% coefficient of pressure vectors
% 
% The function outputs to row vectors. The first is the coefficients of
% drag and the other is the coefficients of lift for every angle of attack

%% Create constants
c = 3.5/39.37;
C_n = 0;
C_a = 0;
aoa = aoa * (pi/180);
% Get the delta x and delta y vector from the file
for j = 1: length(x_pos)-1
    delta_x(j,:) = x_pos(j+1) - x_pos(j);
    delta_y(j,:) = y_pos(j+1) - y_pos(j);
end
delta_x = delta_x/39.37;
delta_y = delta_y/39.37;
%% Calculate the lift and drag coefficients
% Loop through this procedure for every angle of attack
for i = 1: length(C_p)-1
    % Compute the normal and axial coefficients
    C_n = C_n + 0.5*(C_p(:,i) + C_p(:,i+1))*(delta_x(i)/c);
    C_a = C_a + 0.5*(C_p(:,i) + C_p(:,i+1))*(delta_y(i)/c);
end
% Compute the drag coefficient and append it to matrix;
C_n = -C_n;
C_l= C_n*cos(aoa) - C_a*sin(aoa);
C_d= C_n*sin(aoa) + C_a*cos(aoa);
