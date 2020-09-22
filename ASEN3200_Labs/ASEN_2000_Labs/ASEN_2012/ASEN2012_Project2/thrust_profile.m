function thrust = thrust_profile(y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% thrust_profile function
% Created by: Johnathan Tucker
% Author ID #: 108602275
%
% Purpose: to calculate the thrust profile of a water bottle rocket
%
% Inputs: the solution matrix from ode45 calculated in main script
%
% Outputs: the thrust profile for a water bottle rocket
%
% Assumptions: 
%   - Static stability
%   - No wind
%   - Forces due to rail are ignored
%
% Date Created: 11/29/2018
% Date Modified: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create all needed constants
C_discharge = 0.8;
Vol_bottle = 0.002;
P_atm = 83426.56;
gamma = 1.4;
d_throat = .021;
R = 287;
P_initial = 344738 + P_atm;
Vol_water_initial = 0.001;
T_initial = 300;
A_throat = pi*(d_throat/2)^2;
Vol_air_initial = Vol_bottle - Vol_water_initial;
m_air_initial = (P_initial*Vol_air_initial)/(R*T_initial);
%% Start with phase 1
% First we need to find the index that phase one ends
index_change = find(y(:,7) > Vol_bottle, 1);
% Now we can calculate the pressure for these indices
press = ((Vol_air_initial./y(1:index_change,7)).^gamma).*P_initial;
% And the thrust for phase 1
thrust_1 = 2*C_discharge*(press-P_atm)*A_throat;
%% Now phase 2
% First we need the index where it goes from phase 2 to 3
index_change_2 = find(y(index_change:end,2) < m_air_initial, 1, 'last');
% Now we create the press end variable
P_end_2 = P_initial*((Vol_air_initial/Vol_bottle)^gamma);
% Create the mass of air vector from phase 1 index to phase 3 start index
m_air_2 = y(index_change:index_change_2,2);
% calculate the air pressure for phase 2
press_2 = ((m_air_2./m_air_initial).^gamma).*P_end_2;
% calculate the density for phase 2
rho_2 = m_air_2./Vol_bottle;
% calculate the temperature for phase 2
temp_2 = press_2./(rho_2*R);
% calculate the critical pressure for phase 2
P_crit_2 = press_2.*((2/(gamma+1))^(gamma/(gamma-1)));
% loop through the critical pressure vector
for i = 1:length(P_crit_2)
    % if it's choked
    if P_crit_2(i) > P_atm
        % Calculate the necessary exit variables 
        temp_end_2 = (2/(gamma+1))*temp_2(i);
        v_exh_2 = sqrt(gamma*R*temp_end_2);
        rho_exit_2 = P_crit_2(i)/(R*temp_end_2);
        p_exit_2 = P_crit_2(i);
    % if it isn't chokes
    else
        % calculate the mach
        mach_e_2 = sqrt((((press_2(i)/P_atm)^((gamma-1)/gamma))-1)*(2/(gamma-1)));
        % and then the same end variables as choked
        temp_end_2 = temp_2(i)/(1+((gamma-1)/2)*mach_e_2^2);
        v_exh_2 = mach_e_2*sqrt(gamma*R*temp_end_2);
        rho_exit_2 = P_atm/(R*temp_end_2);
        p_exit_2 = P_atm;
    end
    % calculate the mass flow rate
    m_dot_2 = C_discharge*rho_exit_2*A_throat*v_exh_2;
    % then the thrust
    f_thrust = (m_dot_2*v_exh_2) + (p_exit_2 - P_atm)*A_throat;
    % check to make sure there aren't any negatives
    if f_thrust < 0
        f_thrust = 0;
    end
    % create the thrust phase 2 vector and append the forces to it
    thrust_2(i,1) = f_thrust;
end
% there is no thrust in phase 3
thrust_3 = zeros(length(y(index_change_2:end,1))-2,1);
% create the output vector
thrust = [thrust_1; thrust_2; thrust_3];
end

