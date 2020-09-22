function dydt = odefunc(t,y,c_forces, c_moments,flag,controlled)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
dydt = zeros(12,1);
% Flag for nonlinear model
if flag
    %% Create constants
    g = -9.81;
    m = 0.068;
    alpha = (2*10^-6);
    zeta = 1*(10^-3);
    
    % Pull data from input vector
    phi = y(1);
    theta = y(2);
    psi = y(3);
    p = y(4);
    q = y(5);
    r = y(6);
    u = y(10);
    v = y(11);
    w = y(12);
    % Calculate aerodynamic forces
    a_forces = (-zeta*norm([u; v; w]))*[u; v; w];
    a_moments = (-alpha*norm([p;q;r]))* [p;q;r];
    % Check to see if feedback control should be implemented
    if controlled
        c_moments(1) = -p*0.004;
        c_moments(2) = -q*0.004;
        c_moments(3) = -r*0.004;
    end
    G_b = a_moments + c_moments;
    % Calculate the Inertia tensor
    I_x = 6.8*(10^-5);
    I_y = 9.2*(10^-5);
    I_z = 1.35*(10^-4);
    I_mat = [I_x,I_y,I_z];
    I_mat = diag(I_mat);
    %% Begin solving for their time derivatives via dynamics
    euler_rate = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);...
                  0, cos(phi), -sin(phi);...
                  0, sin(phi)*sec(theta), cos(phi)*sec(theta)];
    euler_rate = euler_rate * [p; q; r];
    
    % Calculate angular acceleration vector
    ang_accel = I_mat*[p; q; r];
    ang_accel = cross(-[p; q; r],ang_accel) + G_b;
    ang_accel = inv(I_mat)*ang_accel;
    
    % Get the gravity in body coordinates
    g_b = [-sin(theta); sin(phi)*cos(theta); cos(phi)*cos(theta)].*g;
    
    % Calculate accceleration vector
    accel = g_b + a_forces.*(1/m) + (c_forces).*(1/m) - ...
        cross([p;q;r],[y(10);y(11);y(12)]);

    % Apply transformation matrix to get the velocity in inertial
    % coordinates
    L_BE = angle2dcm(psi,theta,phi);
    V_E = L_BE' * [y(10); y(11); y(12)];
    
    %% Create the output vector
    dydt(1) = euler_rate(1);
    dydt(2) = euler_rate(2);
    dydt(3) = euler_rate(3);
    dydt(4) = ang_accel(1);
    dydt(5) = ang_accel(2);
    dydt(6) = ang_accel(3);
    dydt(7) = V_E(1);
    dydt(8) = V_E(2);
    dydt(9) = V_E(3);
    dydt(10) = accel(1);
    dydt(11) = accel(2);
    dydt(12) = accel(3);
else
    % Create constants
    m = 0.068;
    g = 9.81;
    k = .0024;
    r = 0.060;
    I_x = 6.8*(10^-5);
    I_y = 9.25*(10^-5);
    I_z = 1.35*(10^-5);
    
    % Create deviation constants
    dev_phi = y(1);
    dev_theta = y(2);
    dev_psi = y(3);
    dev_p = y(4);
    dev_q = y(5);
    dev_r = y(6);
    dev_u_E = y(10);
    dev_v_E = y(11);
    dev_w_E = y(12);
    
    % Transform the velocity's from body to inertial coordinates
    L_BE = angle2dcm(dev_psi,dev_theta,dev_phi);
    dev_V_E = L_BE' * [y(10); y(11); y(12)];
    
    % Create output vector
    dydt(1) = dev_p;
    dydt(2) = dev_q;
    dydt(3) = dev_r;
    dydt(4) = (r/(sqrt(2)*I_x))*(0);
    dydt(5) = (r/(sqrt(2)*I_y))*(0);
    dydt(6) = (k/(I_z))*(0);
    dydt(7) = dev_V_E(1);
    dydt(8) = dev_V_E(2);
    dydt(9) = dev_V_E(3);
    dydt(10) = -g*dev_theta;
    dydt(11) = g*dev_phi;
    dydt(12) = (1/m)*(0);
end
end
