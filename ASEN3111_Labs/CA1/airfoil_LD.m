function [L,D] = airfoil_LD(const, num_panels, Cp_upper, Cp_lower)
%AIRFOIL_LD Output lift and drag of the airfoil using trapezoidal rule
%
% Author: Johnathan Tucker
% Collaborators: N/A
% This function takes in the struct of constants, number of panels to use,
% and the upper/lower Cp data. It outputs the lift and drag for the airfoil
% after using the trapezoidal rule to find the normal and axial forces, 
% which are then converted to L and D.
%
% Last Revised: 2/2/2020
%% Begin the trapezoidal rule process for the airfoil
% Create thickness as a fraction of chord
t = 12/100;
% Create vector for percentage of chord length
x = linspace(0,const.chord,num_panels);
% Create a vector of half thicknesses
y = (t/0.2)*const.chord*(0.2969.*sqrt(x./const.chord) - ...
    0.1260.*(x./const.chord) - 0.3516.*(x./const.chord).^2 +...
    0.2843.*(x./const.chord).^3 - 0.1036.*(x./const.chord).^4);

% Create summation variables
N_sum = 0;
A_sum = 0;
% Loop through the number of panels given minus one so that I can fully
% increment through the x vector
for i = 1:(num_panels -1)
    % Calculate the upper and lower pressures for the even and odd increments
    p_u_even = fnval(Cp_upper,x(i+1)/const.chord)*const.q_inf + const.p_inf;
    p_u_odd = fnval(Cp_upper,x(i)/const.chord)*const.q_inf + const.p_inf;
    p_l_even = fnval(Cp_lower,x(i+1)/const.chord)*const.q_inf + const.p_inf;
    p_l_odd = fnval(Cp_lower,x(i)/const.chord)*const.q_inf + const.p_inf;
    
    % Calculate the normal force 
    N  = 0.5*(-(p_u_even + p_u_odd) + (p_l_even + p_l_odd))*(x(i+1) - x(i));
    % Add the normal force to the ongoing summation
    N_sum = N_sum + N;
    % Calculate the axial force 
    A = 0.5*((p_u_even + p_u_odd) + (p_l_even + p_l_odd))*(y(i+1) - y(i));
    % Add the axial force to the ongoing summation
    A_sum = A_sum + A;
end

% Convert from N and A to L and D
L = N_sum*cosd(const.aoa) - A_sum*sind(const.aoa);
D = A_sum*cosd(const.aoa) + N_sum*sind(const.aoa);
end

