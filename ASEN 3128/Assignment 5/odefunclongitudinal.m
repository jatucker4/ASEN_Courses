function dydt = odefunclongitudinal(t,y,c_forces, c_moments,flag,controlled,dev_v_ref,...
    dev_u_ref)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
dydt = zeros(12,1);
% Flag for nonlinear model
% Set gain constants
k1 = 0.001276;
k2 = 0.00232;
k3 = 0.046984899496650;

k4 = 0.001584;
k5 = 0.00288;
k6 = -0.046984899496650;
if flag
    %% Create constants
    g = 9.81;% gravity [m/s^2]
    m = 0.068;% mass [kg]
    alpha = (2*10^-6);% aero moment const [Nm/(rad/s)^2]
    zeta = 1*(10^-3);% aero force const [Nm/(m/s)^2]
    
    % Pull data from input vector
    phi = y(1); % bank angle [rad]
    theta = y(2);% elevation angle [rad]
    psi = y(3);% azimuth angle [rad]
    p = y(4);% roll rate [rad/s]
    q = y(5);% pitch rate [rad/s]
    r = y(6);% yaw rate [rad/s]
    u = y(10);% rel wind x [m/s]
    v = y(11);% rel wind y [m/s]
    w = y(12);% rel wind z [m/s]
    % Calculate aerodynamic forces
    a_forces = (-zeta*norm([u; v; w]))*[u; v; w];% aerodynamic force [N]
    a_moments = (-alpha*norm([p;q;r]))* [p;q;r]; % aerodynamic moment [Nm]
    % Check to see if feedback control should be implemented
    if controlled
          % If it's controlled use the gains calculated in lab
          c_moments(1) = (-0.001276*p-0.00232*phi); % control moment [Nm]
          c_moments(2) = (-0.001584*q-0.00288*theta); % control moment [Nm]
          c_moments(3) = c_moments(3)-r*0.004; % control moment [Nm]           
    end
    G_b = a_moments + c_moments; % sum of moments [Nm]
    % Calculate the Inertia tensor
    I_x = 5.8*(10^-5); % moment of inertia [kg m^2]
    I_y = 7.2*(10^-5); % moment of inertia [kg m^2]
    I_z = 1.0*(10^-4); % moment of inertia [kg m^2]
    I_mat = [I_x,I_y,I_z]; % moment of inertia [kg m^2]
    I_mat = diag(I_mat); % moment of inertia [kg m^2]
    %% Begin solving for their time derivatives via dynamics
    euler_rate = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);...
                  0, cos(phi), -sin(phi);...
                  0, sin(phi)*sec(theta), cos(phi)*sec(theta)];
    euler_rate = euler_rate * [p; q; r]; % change in euler angles [rad/s]
    
    % Calculate angular acceleration vector
    ang_accel = I_mat*[p; q; r]; % angular acceleration [rad/s^2]
    ang_accel = cross(-[p; q; r],ang_accel) + G_b; % angular acceleration [rad/s^2]
    ang_accel = inv(I_mat)*ang_accel; % angular acceleration [rad/s^2]

    % Get the gravity in body coordinates
    g_b = [-sin(theta); sin(phi)*cos(theta); cos(phi)*cos(theta)].*g;
    
    % Calculate accceleration vector
    accel = g_b + a_forces.*(1/m) + (c_forces).*(1/m) - ...
        cross([p;q;r],[y(10);y(11);y(12)]);% acceleration [m/s^2]

    % Apply transformation matrix to get the velocity in inertial
    % coordinates
    L_BE = angle2dcm(psi,theta,phi);% transformation matrix [rad]
    V_E = L_BE' * [y(10); y(11); y(12)];% inertial velocity [m/s]
    
    %% Create the output vector
    dydt(1) = euler_rate(1); % change in euler angles [rad/s]
    dydt(2) = euler_rate(2);% change in euler angles [rad/s]
    dydt(3) = euler_rate(3);% change in euler angles [rad/s]
    dydt(4) = ang_accel(1);% angular acceleration [rad/s^2]
    dydt(5) = ang_accel(2);% angular acceleration [rad/s^2]
    dydt(6) = ang_accel(3);% angular acceleration [rad/s^2]
    dydt(7) = V_E(1);% inertial velocity [m/s]
    dydt(8) = V_E(2);% inertial velocity [m/s]
    dydt(9) = V_E(3);% inertial velocity [m/s]
    dydt(10) = accel(1);% acceleration [m/s^2]
    dydt(11) = accel(2);% acceleration [m/s^2]
    dydt(12) = accel(3);% acceleration [m/s^2]
else
    % Create constants
    m = 0.068; % mass [kg]
    g = 9.81; % gravity [m/s^2]
    k = .0024; % moment constant [m]
    r = 0.060; % distance to thrust vector [m]
    I_x = 6.8*(10^-5); % moment of inertia [kg m^2]
    I_y = 9.25*(10^-5); % moment of inertia [kg m^2]
    I_z = 1.35*(10^-5); % moment of inertia [kg m^2]
    
    % Create deviation constants
    dev_phi = y(1);% deviation in bank angle [rad]
    dev_theta = y(2);% deviation in elevation angle [rad]
    dev_psi = y(3);% deviation in azimuth angle [rad]
    dev_p = y(4);% deviation in roll rate [rad/s]
    dev_q = y(5);% deviation in pitch rate[rad/s]
    dev_r = y(6);% deviation in yaw rate[rad/s]
    dev_u_E = y(10);% deviation in inertial x vel[m/s]
    dev_v_E = y(11);% deviation in inertial y vel[m/s]
    dev_w_E = y(12);% deviation in inertial z vel[m/s]
    
    % Transform the velocity's from body to inertial coordinates
    L_BE = angle2dcm(dev_psi,dev_theta,dev_phi);
    dev_V_E = L_BE' * [y(10); y(11); y(12)];
    
    % Create output vector
    dydt(1) = dev_p;% deviation in roll rate
    dydt(2) = dev_q;% deviation in pitch rate
    dydt(3) = dev_r;% deviation in yaw rate
    % Check to see if the system is controleld
    if controlled
        dydt(4) = 0*((-k1*dev_p-k2*dev_phi - k2*k3*dev_v_E)/(I_x) +...
            (k2*k3*dev_v_ref/I_x));
        dydt(5) = (-k4*dev_q-k5*dev_theta - k5*k6*dev_u_E)/(I_y) +...
            (k5*k6*dev_u_ref/I_y);
        dydt(6) = 0*-0.004*dev_r/I_z;
    else
        dydt(4) = (r/(sqrt(2)*I_x))*(0);% change in roll rate [rad/s^2]
        dydt(5) = (r/(sqrt(2)*I_y))*(0);% change in pitch rate [rad/s^2]
        dydt(6) = (k/(I_z))*(0);% change in yaw rate [rad/s^2]
    end
    dydt(7) = dev_V_E(1);% deviation in inertial vel [m/s]
    dydt(8) = dev_V_E(2);% deviation in inertial vel [m/s]
    dydt(9) = dev_V_E(3);% deviation in inertial vel [m/s]
    dydt(10) = -g*dev_theta;% deviation in vel [m/s]
    dydt(11) = g*dev_phi;% deviation in vel [m/s]
    dydt(12) = (1/m)*(0);% deviation in vel [m/s]
end
end