function [pos_final,arc_length, time,z_t_vec,...
    Gs_felt] = Lab1_Loop(loop_radius,t_start, pos_start, coaster_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by: Johnathan Tuccker in collaboration with Cameron Turman
% Purpose: Create loop piece for the coaster track and output
%          necessary information for other pieces and calculations
% Inputs: loop_radius: Radius of the banked half circle in meters
%         pos_start: The x,y,z coordinate of the beginning of the banked
%                    turn in meters
%         coaster_plot: The figure name to plot the banked portion
%         t_start: Start time of the loop in seconds
% Outputs: pos_final: Final x,y,z coordinates at the end point of the
%                     banked turn in meters
%          arc_length: Length of the banked turn in meters
%          time: Time it took to complete the banked turn in seconds
%          z_t_vec: The change in height over time of the banked turn in
%                   meters
%          GS_felt: Vector of G's for loop section. Unitless
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initial Setup %%
    t_end = 2*pi;

    MIN_T_INTERVAL = t_start;
    v0 = sqrt(2* 9.81*40);
    T = 2*pi*loop_radius/v0;
    omega = 2*pi/T;
    d_off = 1;
    phi = pi/2;
    m = 1;
    % Defining start location of Loop
    tspan =  linspace(0,round(T)+1,112);
    tspan = tspan(tspan <= T);
    x_start = pos_start(1);
    y_start = pos_start(2);
    z_start = pos_start(3)+20;
    
    %% make parametrization of transition %%
    x_t = x_start - loop_radius*cos(omega*tspan+ phi);
    y_t = y_start + (1/T)*d_off*tspan;
    z_t = z_start - loop_radius*sin(omega*tspan + phi);
    z_t_vec = z_t;
    %% Get the final position values for output
    pos_final(1) = x_t(end);
    pos_final(2) = y_t(end);
    pos_final(3) = z_t(end);
    %% Plot
    figure(coaster_plot)
    plot3(x_t, y_t, z_t,'LineWidth',2)
    hold on
    scatter3(pos_final(1),pos_final(2),pos_final(3),'filled')
    %% Calculate arc length
    syms t
    x_t = x_start - loop_radius*cos(omega*t+ phi);
    y_t = y_start + (1/T)*d_off*t;
    z_t = z_start - loop_radius*sin(omega*t + phi);
    s_t = matlabFunction(sqrt((diff(x_t)^2) + (diff(y_t)^2) +...
        (diff(z_t)^2)));
    arc_length = integral(s_t ,0,tspan(end));
    
    %% Find the G's
    % Find Angle of Track in z-direction with respect to horizontal axis
    dxdt = omega*loop_radius*sin(omega *tspan + phi);
    dydt = d_off/T;
    dzdt = -omega*loop_radius*cos(omega *tspan + phi);

    dzdxdy = dzdt./(sqrt(dxdt.^2+dydt.^2)); % slope of line in 3d space

    theta = atan(dzdxdy);

    % Balance Forces to Find Normal Force
    v = sqrt(2*9.81*(125-z_t_vec));
    N = (v.^2/loop_radius)*m + cos(theta)*m*9.81;

    % Calculate G's along
    Gs_felt = N/9.81;
    frontgs = -sin(theta);
    time = 0;
end

