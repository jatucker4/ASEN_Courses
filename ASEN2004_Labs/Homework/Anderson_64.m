% This script is for problem 6.4 in Anderson's introduction to flight
% Create Variables
rho = 0.002377;
S = 181;
k = 1/(pi*0.91*6.2);
P_A_max = .83*345;
V_vec = 1:.1:500;

syms V
q = 0.5*rho*V^2;
C_L = 3000/(q*S);
C_D = 0.027 + k*C_L^2;
T_R = 3000/(C_L/C_D);
P_R = (T_R*V)/550;

P_R_vec = double(subs(P_R,V,V_vec));
P_A_max_vec = P_A_max*ones(1,length(V_vec));
[M_1,I] = min(P_R_vec);

figure(1)
plot(V_vec,P_R_vec)
hold on
plot(V_vec,P_A_max_vec)
title('$Power\:Required\:versus\:Velocity\:and\:Power\:Allowed$',...
    'Interpreter','latex')
xlabel('$Velocity\:(m/s)$','Interpreter','latex')
ylabel('$Power\:(hp)$','Interpreter','latex')
legend('Power Required', 'Power Allowed')

V_max_sea_level = double(solve(P_R == P_A_max));
V_max_sea_level = V_max_sea_level(2);
fprintf("The velocity max is %f at sea level\n",V_max_sea_level);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now for altitude
rho_alt = 0.001648;

V_alt = (rho/rho_alt)^0.5 * V;
P_R_alt = (rho/rho_alt)^0.5 * P_R;
P_A_alt = (rho/rho_alt)^-1 * P_A_max;

V_alt_vec = double(subs(V_alt,V,V_vec));
P_R_alt_vec = double(subs(P_R_alt,V,V_vec));
P_A_alt_vec = P_A_alt*ones(1,length(V_vec));
[M,I] = min(P_R_alt_vec);


figure(2)
plot(V_alt_vec,P_R_alt_vec)
hold on
plot(V_alt_vec,P_A_alt_vec)

title('$Power\:Required\:versus\:Velocity\:and\:Power\:Allowed\:at\:Altitude$',...
    'Interpreter','latex')
xlabel('$Velocity\:(m/s)$','Interpreter','latex')
ylabel('$Power\:(hp)$','Interpreter','latex')
legend('Power Required', 'Power Allowed')

tol = 0.1; % change tolerance as required.
index = abs(P_R_alt_vec-P_A_alt_vec) < tol;
indexNumeric = find(index);
V_max_alt = V_alt_vec(indexNumeric(2));

fprintf("The velocity max is %f at altitude\n",V_max_alt);
