function [r_c, v_c] = oe2cart(oe_c, dt , mu);
%OE2CART Summary of this function goes here
%   Detailed explanation goes here
a_c = oe_c(1);
e_c = oe_c(2);
i_c = oe_c(3);
Ohm_c = oe_c(4);
w_c = oe_c(5);
M0_c = oe_c(6);

syms t

eqn = M0_c == 2*atan(sqrt((1-e_c)/(1+e_c))*tan(t/2)) - ...
    (e_c*sqrt(1-e_c^2)*sin(t))/(1+e_c*cos(t));

t_anomaly = solve(eqn,t);
t_anomaly = t_anomaly*180/pi;

x = a_c*(e_c + cos(t_anomaly))/(1+e_c*cos(t_anomaly));
y = a_c*(1 - e_c^2)*sin(t_anomaly)/(1 + e_c*cos(t_anomaly));

r_c = sqrt(x^2 + y^2);

v_c = sqrt(2*((-mu/(2*a_c)) + mu/r_c));

end

