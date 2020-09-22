%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear;
close all;
clc;
%% Create variables
[T, a, P, rho] = atmoscoesa(1624, 'None');
% Chord Length
c = .064;
% Wing span
b = .70;
% Planform Area
S = c*b;
% Aspect ratio
AR = b^2 / S;
% Create a velocity vector for required velocity range
V = 3:0.5:7;
e = 0.8;
k = 1/(pi * e * AR);
%% Calculate Swet
wing_wet = 2*1.25*0.7*.064;
wing_wet_total = 2*wing_wet;
tail_wet_total = 2*1.25*.15*.05;
fuselage_length = 0.30;
fuselage_radius = .025;
fuselage_wet = 2*pi*fuselage_radius*fuselage_length + 2*pi*fuselage_radius^2;
rudder_wet = .08*.05;
S_wet = fuselage_wet + wing_wet_total + tail_wet_total + rudder_wet;
%% Calculate weight
rho_balsa = 160;
rho_foam = 40;
balsa_thickness = .0127;
rho_super_glue = 1.1;
rho_tissue_paper = 24; % note this is g/m^2
W_fuselage = pi * fuselage_radius^2 * fuselage_length *rho_foam;
W_horz_tail = rudder_wet * balsa_thickness * rho_foam;
W_vert_tail = tail_wet_total * balsa_thickness * rho_foam;
W_wing = S * 0.005564 * rho_balsa + (2*1.25*S*rho_tissue_paper)/1000;
W_ballast = 0.003;

W_total = W_fuselage + W_horz_tail + W_vert_tail + W_wing + W_ballast;

%% Calculate CG
x_cg_fuselage = .0475;
x_cg_wing = .064 * .40;
x_cg_tail = .08 + .095;
x_cg_horz_tail = .05 + .095 + .025;
x_cg_ballast = -.175;

% x_cg_fuselage = .0625;
% x_cg_wing = .064 * .30;
% x_cg_tail = .025 + .125;
% x_cg_horz_tail = .05 + .125;
% x_cg_ballast = -.175;


x_cg_tot = (W_fuselage*x_cg_fuselage + W_horz_tail*x_cg_horz_tail +...
    W_wing*x_cg_wing + W_vert_tail*x_cg_tail + W_ballast*x_cg_ballast)/...
    W_total;

wing_loading = 9.81 * W_total/S;
% Calculate desired L/Dmax
L_D_desired = 50/8;

%% Calculate Cd0
Re_temp = rho * 6 * 0.35 / (1.74 * 10^-5);
C_fe = 0.074 / Re_temp^0.2;
Cd0 = C_fe * S_wet / S;


%% Calculate CL
CL = L_D_desired * 2 * Cd0;

% Calculate required wing loading
Wing_loading_req = .3022 * 0.5 * rho * V.^2 ;

plot(V,Wing_loading_req);
hold on 
plot(V,ones(1,length(Wing_loading_req))*7.189)
title('$Wing\:Loading\:vs\:Velocity$','Interpreter','latex')
xlabel('$Velocity\:(m/s)$','Interpreter','latex')
ylabel('$Wing\:Loading\:(Pa)$','Interpreter','latex')
legend('Wing Loading Calculated','Wing Loading Required')
%% Calculate VH and VV
Vh = (.15*.05) * (x_cg_horz_tail - x_cg_tot ) / (S*c);
Vv = (.08*.05) * (x_cg_tail - x_cg_tot ) / (b*c);

% Glide range velocity
V_req = sqrt(2*W_total*9.81 / (rho*S*sqrt(Cd0/k)));