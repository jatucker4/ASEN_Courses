function [pos_final,arc_length, time, ...
    z_t_vec] = Lab1_Helix(helix_radius,helix_loops,helix_height,...
    t_start, pos_start, coaster_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by: Johnathan Tucker
% Purpose: Create Helix piece for the coaster track and output
%          necessary information for other pieces and calculations
% Inputs: helix_radius: Radius of the helix circle in meters
%         helix_height: Height of the helix in meters
%         t_start: Time the helix started in seconds
%         pos_start: The x,y,z coordinate of the beginning of the banked
%                    turn in meters
%         coaster_plot: The figure name to plot the banked portion
%
% Outputs: pos_final: Final x,y,z coordinates at the end point of the
%                     banked turn in meters
%          arc_length: Length of the banked turn in meters
%          time: Time it took to complete the banked turn in seconds
%          z_t_vec: The change in height over time of the banked turn
%                   in meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create variables
    t_end = 2*pi + t_start;
    t_interval = t_end - t_start;
    x_start = pos_start(1);
    y_start = pos_start(2);
    z_start = pos_start(3);
    % Create parametric functions
    syms t s t_new
    x_t = helix_radius * sin(helix_loops * (t - t_start)) + x_start;
    y_t = helix_radius * cos(helix_loops * (t - t_start)) + (y_start - ...
                                                   helix_radius);
    z_t = -helix_height * (t - t_start) / t_interval + z_start;

    %% Get the final position values for output
    pos_final(1) = double(subs(x_t,t,t_end));
    pos_final(2) = double(subs(y_t,t,t_end));
    pos_final(3) = double(subs(z_t,t,t_end));
    
    %% Calculate arc length of the segment
    s_t = matlabFunction(sqrt((diff(x_t)^2) + ...
        (diff(y_t)^2) + (diff(z_t)^2)));
    arc_length = integral(s_t ,t_start,t_end);
    
    %% Plot Helix
    figure(coaster_plot)
    fplot3(x_t, y_t, z_t, [t_start, t_end], 'LineWidth',2)
    time = 0;
    %% get z_vec for output
    tspan = linspace(t_start,t_end);
    z_t_vec = double(subs(z_t,t,tspan));
end

