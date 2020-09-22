function dx_dt = g_fun(t, y, mu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
theta = y(1);
omega = y(2);

dtheta_dt = omega;
domega_dt = mu*(1-(theta^2))*dtheta_dt - theta;

dx_dt = [dtheta_dt; domega_dt]; 
end

