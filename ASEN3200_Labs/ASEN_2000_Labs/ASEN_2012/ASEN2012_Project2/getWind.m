function [v_wind_surface,v_wind_aloft] = getWind(cardinal_surf_wind,...
    cardinal_aloft_wind)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 % [magnitude of wind (m/s), theta (rad), vertical wind (m/s)]

    %%% Surface Wind %%%
    phi = 204 * (pi/180);
    wind_speed = cardinal_surf_wind(1);
    theta = cardinal_surf_wind(2);
    vert_wind_speed = cardinal_surf_wind(3);
    
    % make wind vector in cardinal coord. frame
    v_wind_cardinal = [-wind_speed * cos(theta);
                       wind_speed * sin(theta);
                       vert_wind_speed];
    
    % Build rot. matrix, A, to transform - A: cardinal -> launch pad 
    A = [cos(phi)  -sin(phi) 0;
         sin(phi)   cos(phi) 0;
         0          0        1];
    
    % Transform v_wind_cardinal to v_wind_launch_pad:
    v_wind_surf_launch_pad = A * v_wind_cardinal;



    %%% Aloft Wind %%%
    wind_speed = cardinal_aloft_wind(1);
    theta = cardinal_aloft_wind(2);
    vert_wind_speed = cardinal_aloft_wind(3);
    
    % make wind vector in cardinal coord. frame
    v_wind_cardinal = [-wind_speed * cos(theta);
                       wind_speed * sin(theta);
                       vert_wind_speed];
    
    % Build rot. matrix, A, to transform - A: cardinal -> launch pad 
    A = [cos(phi) -sin(phi) 0;
         sin(phi)  cos(phi) 0;
         0         0        1];
    
    % Transform v_wind_cardinal to v_wind_launch_pad:
    v_wind_aloft_launch_pad = A * v_wind_cardinal;



    % Assigning Global Vars
    % [perp to launch dir, || to ground, launch dir, up]
    v_wind_surface = v_wind_surf_launch_pad; % [m/s, m/s, m/s]
    v_wind_aloft = v_wind_aloft_launch_pad; % [m/s, m/s, m/s]
end

