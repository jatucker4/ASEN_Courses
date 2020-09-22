function dy = isp_function(~,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is defined to be called by ode45. It uses the "Isp model"
% to model the rockets flight trajectory. It's important to note that the
% Isp model assumes that all of the fuel is expended instantaneously
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create constants
g = 9.81;
C_discharge = 0.8;
rho_air_amb = 0.961;
Vol_bottle = 0.002;
P_atm = 83426.56;
gamma = 1.4;
rho_water = 1000;
d_throat = .021;
d_bottle = .105;
R = 287;
m_bottle = 0.15;
C_drag = 0.5;
P_initial = 344738 + P_atm;
Vol_water_initial = 0.001;
T_initial = 300;
rail_length = 0.4;
A_throat = pi*(d_throat/2)^2;
A_bottle = pi*(d_bottle/2)^2;
Vol_air_initial = Vol_bottle - Vol_water_initial;
m_air_initial = (P_initial*Vol_air_initial)/(R*T_initial);
wind_aloft = 0;
wind_surface =0;

%% Get initial values
velocity_x = y(1);
m_air = y(2);
m_water = y(3);
theta = y(4);
x_pos = y(5);
z_pos = y(6);
Vol_air = y(7);

% Get wind effects
v_rel = ((wind_aloft - wind_surface) / 8) * z_pos;
% Must create a P_air in the case that it isn't phase 1
if ~(Vol_air < Vol_bottle)
    % Get the end Pressure
    P_end = P_initial * (Vol_air_initial/Vol_bottle)^gamma;
        
    % Get the Pressure for the second phase
    P_air = P_end * (m_air/m_air_initial)^(gamma);
end

if Vol_air > Vol_bottle
    Vol_air = Vol_bottle;
end
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
    F = 2*C_discharge*A_throat*(P_air - P_atm);
    
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
    F = -dm_air_dt*V_exit + (P_e - P_atm)*A_throat;
    
%% Now the final phase
else
    % Thrust is zero and there is no change in mass or volume
    F = 0;
    dm_water_dt = 0;
    dm_air_dt = 0;
    dV_dt = 0;
end
%% Common equations between all states
% Get drag
Drag = (rho_air_amb/2)*(velocity^2)*C_drag*A_bottle;

% Make sure mass of water doesn't go negative
if m_water < 0
    m_water = 0;
end

% Calculate mass of the rocket
m_rocket = m_bottle + m_water + m_air;

% Calculate change in velocity
dvelocity_dt = (F - Drag - (m_rocket*g*sin(theta)))/m_rocket;

% Get the change in the x and z positions
dx_dt = cos(theta) * velocity;
dz_dt = sin(theta) * velocity;

% Get the change in angle taking into account the stand
if rail_length * cos(theta) < z_pos
    
    dtheta_dt = (-g*cos(theta))/velocity;

else
    dtheta_dt = 0;
end
% Account for wind effects after research I found that surface boundary
% layer exists until about 1 to 2 km


%% Now output everything to ode45
dy = zeros(7,1);
dy(1) = dvelocity_dt;
dy(2) = dm_air_dt;
dy(3) = dm_water_dt;
dy(4) = dtheta_dt;
dy(5) = dx_dt;
dy(6) = dz_dt;
dy(7) = dV_dt;


end

