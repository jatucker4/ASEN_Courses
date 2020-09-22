function [pos_final, arc_length, time, z_t_vec, ...
    k_avg] = transition(pos_start, coaster_plot, flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by: Johnathan Tucker
% Purpose: Create transition pieces for the coaster track and output
%          necessary information for other pieces and calculations
% Inputs: flag: Determines which transition to use
%         pos_start: The x,y,z coordinate of the beginning of the banked
%                    turn in meters
%         coaster_plot: The figure name to plot the banked portion
%
% Outputs: pos_final: Final x,y,z coordinates at the end point of the
%                     banked turn in meters
%          arc_length: Length of the banked turn in meters
%          time: Time it took to complete the banked turn in seconds
%          z_t_vec: The change in height over time of the banked turn in
%                   meters
%          k_avg: average radius of curvature in meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_start = 0;
    t_end = 1.881473;
   
    % Defining start location of Helix 
    x_start = pos_start(1);
    y_start = pos_start(2);
    z_start = pos_start(3);
    
    
    %% make parametrization of transition %%
    if flag == 0
        syms t
        x_t = x_start + 20*t;
        y_t = y_start;
        z_t = z_start - 86.9*t/15;
        s_t = @(t) sqrt((20)^2 + (-86.9/15)^2);
        arc_length = t_end * s_t(0);
        tspan = linspace(t_start,t_end);
        z_t_vec = double(subs(z_t,t,tspan));
        k_avg = 0;
    else
        syms t
        x_t = (t) + x_start;
        y_t = y_start;
        z_t = (1/250) * (t - (25 * 2))^2 + z_start/2 + 32;
        s_t = matlabFunction(sqrt(1 + (diff(z_t)^2)));
        t_end = 65.81139;
        arc_length = integral(s_t,t_start,t_end);
        tspan = linspace(t_start,t_end);
        z_t_vec = double(subs(z_t,t,tspan));
        r_t = [x_t,y_t,z_t];
        v_t = diff(r_t,t);
        T = v_t./norm(v_t);
        k = diff(T,t)./v_t;
        k_vec = double(subs(k,t,tspan));
        for i = 1: length(k_vec)
            if isnan(k_vec(i))
                k_vec(i) = 0;
            end
        end
        k_avg = mean(k_vec)*-1;
    end
    
    %% Get the final position values for output
    pos_final(1) = double(subs(x_t,t,t_end));
    pos_final(2) = double(subs(y_t,t,t_end));
    pos_final(3) = double(subs(z_t,t,t_end));
  
    %% Plot transitions
    y_t = @(t) y_start;
    figure(coaster_plot)
    fplot3(x_t, y_t, z_t, [t_start, t_end], 'LineWidth', 2)
    time = 0;
end

