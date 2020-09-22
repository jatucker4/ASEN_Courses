function [pos_final,arc_length, time,z_t_vec, ...
    gs_felt] = Lab1_Bank(bank_radius,pos_start, coaster_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by: Johnathan Tucker in collaboration with Justin Troche
% Purpose: Create bank piece for the coaster track and output
%          necessary information for other pieces and calculations
% Inputs: bank_radius: Radius of the banked half circle in meters
%         pos_start: The x,y,z coordinate of the beginning of the banked
%                    turn in meters
%         coaster_plot: The figure name to plot the banked portion
%
% Outputs: pos_final: Final x,y,z coordinates at the end point of the
%                     banked turn in meters
%          arc_length: Length of the banked turn in meters
%          time: Time it took to complete the banked turn
%          z_t_vec: The change in height over time of the banked turn in
%                   meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create necessary variables
t = linspace( -pi/2, pi/2, 100);
g = 9.81;
h0 = 125;
h = 75.0334;
% Velocity function
v = sqrt(2*g*(h0 - h));
% Angle the path is banked at
bank_angle = rad2deg(atan((v^2)/(bank_radius*g)));
% Calculate g's felt
gs_felt = cosd(bank_angle) + ((v^2)/(bank_radius*g))*sind(bank_angle);
% Calculate arc length
arc_length = ((pi*2*bank_radius)/2)+2*bank_radius;

% Initiate x,y,z positions
x_initial = pos_start(1);
y_initial = pos_start(2) + bank_radius;
z_initial = pos_start(3);

% Parametric forumals of banked turn
x = bank_radius*cos(t) + x_initial;
y = bank_radius*sin(t) + y_initial;
z = zeros(1,numel(x)) + z_initial;
z_t_vec = z;
% Get the final position
pos_final = [x(end),y(end),z(end)];
time = 5;
% Plot the coaster
figure(coaster_plot)
plot3(x,y,z, 'LINEWIDTH', 2)
end
