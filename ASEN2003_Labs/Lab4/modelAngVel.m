function [modOmega] = modelAngVel(theta)

M = 11.7; % mass of cylinder (kg)
M0 = 0.7; % mass of trailing supports (kg)
m = 3.4; % mass of extra mass (kg)
R = 0.235; % radius of cylinder (m)
k = 0.203; % radius of gyration of wheel (m)
I = M*k^2; % moment of Inertia
beta = 5.5; % slope of ramp (deg)
r = 0.178; % radius to extra mass (m)
rem = 0.019; % radius of extra mass (m)
g = 9.81; % gravitational acceleration (m/s^2)

omega1 = [1:16];
omega2 = [1:16];
omega3 = [1:16];
omega4 = [1:16];

modOmega = [omega1,omega2,omega3,omega4];

figure;
plot(theta,omega1);
hold on;
plot(theta,omega2);
plot(theta,omega3);
plot(theta,omega4);
xlabel('Angular Position Theta (rad)');
ylabel('Angular Velocity Omega (rad/s)');
title('Angular Velocity vs Angular Position (Model)');
legend('Model 1','Model 2','Model 3','Model 4');
hold off;

end