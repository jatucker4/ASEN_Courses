function [arc_length,pos_final_braking,...
    z_t_vec] = Lab1_Brake(pos_start,coaster_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%         pos_start: The x,y,z coordinate of the beginning of the banked
%                    turn. Meters
%         coaster_plot: The figure name to plot the banked portion
%
% Outputs: pos_final: Final x,y,z coordinates at the end point of the
%                     banked turn. Meters
%          arc_length: Length of the banked turn. Meters
%       
%          z_t_vec: The change in height over time of the banked turn.
%                   Meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create necessary variables
t_start = 0;
t_end = 75;
% Create initial positions
x_initial = pos_start(1);
y_initial = pos_start(2);
z_initial = pos_start(3);
% Create parametric functions
syms t
x_t = x_initial + t;
y_t = y_initial;
z_t = z_initial;
% Find the change in heights over time
tspan = linspace(0,20);
z_t_vec = double(subs(z_t,t,tspan));
% Find the final position
pos_final_braking(1) = double(subs(x_t,t,t_end));
pos_final_braking(2) = double(subs(y_t,t,t_end));
pos_final_braking(3) = double(subs(z_t,t,t_end));
% Calculate arc length
arc_length = 75;

y_t = @(t) y_initial;
z_t = @(t) z_initial;
% Plot the braking segment
figure(coaster_plot)
fplot3(x_t, y_t, z_t, [t_start, t_end], '-.', 'LINEWIDTH',2);
title('Roller Coaster Path', 'FontSize',16);
xlabel('x-location(m)','FontSize',16);
ylabel('y-location(m)','FontSize',16);
zlabel('height(m)','FontSize',16);
end

