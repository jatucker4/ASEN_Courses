function dy = Project_2_diffeqs(~,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Project_2_diffeqs function
% Created by: Johnathan Tucker
% Author ID: 108602275
%
% Purpose: This function serves as a "processing function" for ode45 where
% depending on the outputs of ode45 different conditionals, that represent
% the phases of the water bottle rocket flight, will activate and output
% differential equations for ode45 to process.
% 
% Inputs: Initial conditions for:
%           -velocity
%           -m_air
%           -m_water
%           -theta
%           -x_position
%           -y_position
%           -volume
%
% Outputs: Differential equations for:
%           -velocity
%           -m_air
%           -m_water
%           -theta
%           -x_position
%           -y_position
%           -volume
%
% Assumptions: 
%   - Static stability
%   - No wind
%   - Forces due to rail are ignored
%
% Date Created: 11/23/2018
% Date Modified: 11/24/2018
%                11/26/2018
%                11/28/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global m_bottle_temp C_drag_temp theta_temp Vol_water_initial P_initial wind_surface wind_aloft
%% Create constants
g = [0;0;9.81];
C_discharge = 0.8;
rho_air_amb = 0.961;
rho_air_amb = 1.184;
Vol_bottle = 0.002;
P_atm = 83054.2464;
gamma = 1.4;
rho_water = 1000;
d_throat = .021;
d_bottle = .105;
R = 287;
% theta_temp = 45 * (pi/180);
% m_bottle_temp = 0.117;
% C_drag_temp = 0.33;
% P_initial = 275790 + P_atm;
% Vol_water_initial = 0.00962;
T_initial = 300;
rail_length = 0.5;
A_throat = pi*(d_throat/2)^2;
A_bottle = pi*(d_bottle/2)^2;
Vol_air_initial = Vol_bottle - Vol_water_initial;
m_air_initial = (P_initial*Vol_air_initial)/(R*T_initial);

%% Correct wind for direction
% [wind_aloft,wind_surface] = getWind([1,deg2rad(-22.5),0],[3,deg2rad(22.5),0]);
%% Get initial values
x_pos = y(1);
y_pos = y(2);
z_pos = y(3);
vel_x = y(4);
vel_y = y(5);
vel_z = y(6);
m_air = y(7);
m_water = y(8);
Vol_air = y(9);

% Get wind effects
if z_pos < 8
    vrel = [vel_x,vel_y,vel_z]' - wind_surface;
    vrel_mag = norm(vrel);
else
    vrel = [vel_x,vel_y,vel_z]' - wind_aloft;
    vrel_mag = norm(vrel);
end

% Get the change in angle taking into account the stand
if norm([x_pos,y_pos,z_pos]) < rail_length
    heading = [cos(theta_temp); 0 ; sin(theta_temp)]./...
        norm([cos(theta_temp); 0; sin(theta_temp)]);
else
    heading = (vrel * 1/vrel_mag);
%     heading = vrel;
end

if Vol_air > Vol_bottle
    Vol_air = Vol_bottle;
end
% Must create a P_air in the case that it isn't phase 1
if ~(Vol_air < Vol_bottle)
    % Get the end Pressure
    P_end = P_initial * (Vol_air_initial/Vol_bottle)^gamma;
    % Get the Pressure for the second phase
    P_air = P_end * (m_air/m_air_initial)^(gamma);
end
% Get drag
Drag = ((rho_air_amb/2)*(vrel_mag^2)*C_drag_temp*A_bottle).*heading;

%% Start with Phase 1
if Vol_air < Vol_bottle
    % Get the pressure of the air
    P_air = P_initial*(Vol_air_initial/Vol_air)^gamma;
    % Get the exhaust velocity
    V_exhaust = sqrt((2/rho_water)*(P_air - P_atm));
    
    % There is no mass of air flowing out of the bottle only water
    dm_water_dt = -C_discharge*rho_water*A_throat*V_exhaust;
    dm_air_dt = 0; 
    
    % Get the thrust of the rocket
%     F = (2*C_discharge*A_throat*(P_air - P_atm)).*heading;
    F = (-dm_water_dt*V_exhaust).*heading;
    % Get the rate of change of the volume
    dV_dt = C_discharge*A_throat*V_exhaust;
%% Start phase two
elseif P_air > P_atm
    % Get the end Temperature
    T_end = T_initial * (Vol_air_initial/Vol_bottle)^(gamma - 1);
    
    % Get the critical pressure
    P_critical = P_air*(2/(gamma + 1))^(gamma/(gamma-1));
    
    % Get the density and temperature
    rho_temp = m_air/Vol_bottle;
    Temp = P_air/(R*rho_temp);
    
    if P_critical > P_atm
        % Get the exit temp
        Temp_e = (2/(gamma + 1))*Temp;
        
        % Get the exit velocity
        V_exit = sqrt(gamma*R*Temp_e);
        
        % Get the exit pressure
        P_e = P_critical;
        
        % Get the exit density
        rho_e = P_e/(R*Temp_e);
    else
        % Get the exit pressure
        P_e = P_atm;
        
        % Get the Mach number
        M_e = sqrt((2/(gamma-1))*(((P_air/P_e)^((gamma-1)/gamma))-1));
        
        % Get the exit temp
        Temp_e = Temp/(1+((gamma-1)/2)*(M_e^2));
        
        % Get the exit density
        rho_e = P_e/(R*Temp_e);
        
        % Get the exit velocity
        V_exit = M_e*sqrt(gamma*R*Temp_e);
    end
    
    % Get the change of mass (only air since all water is gone)
    dm_water_dt = 0;
    dm_air_dt = -C_discharge*rho_e*A_throat*V_exit;
    % Volume is no longer changing
    dV_dt = 0;
    % Get the thrust
    F = (-dm_air_dt*V_exit + (P_e - P_atm)*A_throat).*heading;
    
%% Now the final phase
else
    % Thrust is zero and there is no change in mass or volume
    F = [0,0,0];
    dm_water_dt = 0;
    dm_air_dt = 0;
    dV_dt = 0;
end
%% Common equations between all states
% Make sure mass of water doesn't go negative
if m_water < 0
    m_water = 0;
end

% Calculate mass of the rocket
m_rocket = m_bottle_temp + m_water + m_air;

% Calculate change in velocity
dvelocity_dt = (F - Drag - (m_rocket.*g))./m_rocket;
dvelx_dt = dvelocity_dt(1);
dvely_dt = dvelocity_dt(2);
dvelz_dt = dvelocity_dt(3);

% Get the change in the x and z positions
dx_dt = vel_x;
dy_dt = vel_y;
dz_dt = vel_z;
% Account for wind effects after research I found that surface boundary
% layer exists until about 1 to 2 km


%% Now output everything to ode45
dy = zeros(9,1);
dy(1) = dx_dt;
dy(2) = dy_dt;
dy(3) = dz_dt;
dy(4) = dvelx_dt;
dy(5) = dvely_dt;
dy(6) = dvelz_dt;
dy(7) = dm_air_dt;
dy(8) = dm_water_dt;
dy(9) = dV_dt;

end

