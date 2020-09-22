%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2004 Lab 1
% 
% Created by: Johnathan Tucker in collaboration with the purple cobras
%
% Purpose: The driver script for all calculations and necessary functions. 
%          There are no inputs or outputs.
%           
%          Last Updated: 1/30/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear;
close all;
clc;
%% Declare Variables
efficiency = 0.9;
b = .9;
S_wing = 47564.0309263 / (1000^2) * 2;
AR_UAS = b^2 / S_wing;
%% Create 2D data
% alpha_2d = -5:1:15;
% CL_2d = [-0.2446, -0.1465, -.0401, .0568, .1717, .2737, .4058, .5143,...
%     .6167, .7194, .8201, .9193, 1.0129, 1.1027, 1.1844, 1.2533, 1.2865,...
%     1.2763, 1.2329, 1.1635, 1.0951];
% CD_2d = [0.0140, .0091, .0073, .0059, .0049, .0043, .0045, .0050, .0057,...
%     .0066, .0078, .0092, .0112, .0134, .0165, .0201, .0252, .0332,...
%     .0475, .072, .1052];
A = xlsread('airfoildata.xlsx');
alpha_2d = A(:,1)';
CL_2d = A(:,2)';
CD_2d = A(:,3)';
%% Create UAS data
alpha_UAS = -5:1:12;
CL_UAS = [-.32438, -.21503, -.10081, .010503, .12155, .24163, .34336,...
    .45256, .56037, .66625, .76942, .86923, .96386, 1.0441, 1.0745,...
    1.0807, 1.0379, 1.034];
CD_UAS = [.044251, .033783, .028627, .025864, .024643, .025099, .025635,...
    .02766, .030677, .034855, .040403, .04759, .057108, .070132, .090921,...
    .11193, .13254, .15645];
%% CL plots
CL_3d = C_L_Lab1(alpha_2d,CL_2d,efficiency,AR_UAS);
figure(1)
plot(alpha_2d,CL_2d)
% hold on
% plot(alpha_UAS,CL_UAS)
hold on
plot(alpha_2d,CL_3d)

title('$Coefficient\:of\:Lift\:versus\:Angle\:of\:Attack$',...
    'Interpreter','latex')
xlabel('$Angle\:of\:Attack(degrees)$','Interpreter','latex')
ylabel('$Coefficient\:of\:Lift$','Interpreter','latex')
legend('CL\_2D','CL\_3D')


%% 3D Wing Drag Polar
% 2D Wing Drag Polar
CD_3d = CD_2d + ((CL_3d.^2)/(pi*efficiency*AR_UAS));
figure(2)
subplot(1,2,1)
% plot(CL_2d,CD_2d)
% hold on
plot(CL_3d,CD_3d)
title('$Wing\:Drag\:Polar$','Interpreter','latex')
xlabel('$Coefficient\:of\:Lift$','Interpreter','latex')
ylabel('$Coefficient\:of\:Drag$','Interpreter','latex')
legend('3D Wing')
%% Swet estimation
wing_wet = 2*1.25*0.7*.064;
wing_wet_total = 2*wing_wet;
tail_wet_total = 2*1.25*.1*.03;
fuselage_length = 0.30;
fuselage_radius = .025;
fuselage_wet = 2*pi*fuselage_radius*fuselage_length + 2*pi*fuselage_radius^2;
rudder_wet = .05*.05;
S_wet = 0.2931;
%% Whole aircraft drag polar
[M,I] = min(CD_3d);
I = 17;
alpha_min_cd = alpha_2d(I);
alpha_L_0 = alpha_2d(1,32);

x = alpha_2d(1,48:56);
y = CL_3d(1,48:56);
p = polyfit(x,y,1);
a_0 = p(1);
a = a_0/(1+ (57.3*a_0/(pi*efficiency*AR_UAS)));
CL_min_3d = a*(alpha_min_cd - alpha_L_0);

CD_min = 0.004*(S_wet/S_wing);

e_0 = 1.78*(1 - 0.045*AR_UAS^0.68) - 0.64;
k1 = 1/(pi*e_0*AR_UAS);
k2 = -2*k1*CL_min_3d;
CD_0 = CD_min + k1*(CL_min_3d^2);

CD_whole_body = CD_0 + k1*(CL_3d.^2) + k2*CL_3d;

subplot(1,2,2)
plot(CL_3d,CD_whole_body)
% hold on
% plot(CL_UAS,CD_UAS)
title('$Whole\:Body\:Drag\:Polar$','Interpreter','latex')
xlabel('$Coefficient\:of\:Lift$','Interpreter','latex')
ylabel('$Coefficient\:of\:Drag$','Interpreter','latex')
legend('Whole Body Estimate')

figure(3)
plot(CL_3d,CD_whole_body)
hold on
plot(CL_UAS,CD_UAS)
title('$Whole\:Body\:Drag\:Polar$','Interpreter','latex')
xlabel('$Coefficient\:of\:Lift$','Interpreter','latex')
ylabel('$Coefficient\:of\:Drag$','Interpreter','latex')
legend('Whole Body Estimate','UAS')
%% L/D
L_D_aircraft = CL_3d./CD_whole_body;
L_D_CAD = CL_UAS./CD_UAS;

figure(4)
plot(alpha_2d,L_D_aircraft)
% hold on
% plot(alpha_UAS,L_D_CAD)
title('$Lift\:Over\:Drag\:Plot$','Interpreter','latex')
xlabel('$Angle\:of\:Attack$','Interpreter','latex')
ylabel('$L\:/\:D$','Interpreter','latex')
legend('Whole Aircraft')

%% Calculate CD_0 for UAS Tempest
% Find L/D max
[M, I] = max(L_D_aircraft);
I = 61;
% CD_0 = .0449;
L_D_aircraft_max = L_D_aircraft(I);
[M_1, I_1] = max(L_D_CAD);
L_D_CAD_max = L_D_CAD(I_1);
C_L_range_aircraft = sqrt(CD_0/k1);

[M,I] = min(CD_UAS);
alpha_min_cd_UAS = alpha_UAS(I);
alpha_L_0_UAS = alpha_UAS(1,3) + (0 - CL_UAS(1,3))*(alpha_UAS(1,4)-alpha_UAS(1,3))/(CL_UAS(1,4)-CL_UAS(1,3));

x = alpha_UAS(1,3):alpha_UAS(1,10);
y = CL_UAS(1,3:10);
p = polyfit(x,y,1);
a_0_UAS = p(1);
a = a_0_UAS/(1+ (57.3*a_0_UAS/(pi*efficiency*AR_UAS)));
CL_min_UAS = a*(alpha_min_cd_UAS - alpha_L_0_UAS);

CD_min = 0.004*(S_wet/S_wing);

e_0 = 1.78*(1 - 0.045*AR_UAS^0.68) - 0.64;
k1 = 1/(pi*e_0*AR_UAS);
k2 = -2*k1*CL_min_UAS;
CD_0_UAS = CD_min + k1*(CL_min_UAS^2);
C_L_range_UAS = sqrt(CD_0_UAS/k1);

alpha_range_max_aircraft = (C_L_range_aircraft/a)+alpha_L_0;
alpha_range_max_CAD = (C_L_range_UAS/a)+alpha_L_0_UAS;
%% Performance calculations
% Glide range is given by R = h(L/D)
% So in order to maximize glide range we want to maximize L/D because h is
% a constant value
h = 80; % (m)
GTOW = 0.077*9.81;
rho = 1.0451;
v_range_aircraft = sqrt(2*GTOW/(C_L_range_aircraft*rho*S_wing));
v_range_UAS = sqrt(2*GTOW/(C_L_range_UAS*rho*S_wing));
fprintf("L/D max for the whole aircraft is: %f\n",L_D_aircraft_max);
fprintf("L/D max for the UAS is: %f\n",L_D_CAD_max);
fprintf("Velocity that corresponds to L/D max and max range for the whole aircraft is: %f\n",...
    v_range_aircraft);
fprintf("Velocity that corresponds to L/D max and max range for the UAS is: %f\n",...
    v_range_UAS);
fprintf("Angle of attack for max range of the whole aircraft is: %f\n",...
    alpha_range_max_aircraft);
fprintf("Angle of attack for max range of the UAS is: %f\n",...
    alpha_range_max_aircraft);
% Find CL^(3/2)/CD max
% First whole aircraft:
C_L_endurance_aircraft = sqrt(3*CD_0/k1);
C_L_endurance_UAS = sqrt(3*CD_0_UAS/k1);
alpha_endurance_max_aircraft = (C_L_endurance_aircraft/a)+alpha_L_0;
alpha_endurance_max_UAS = (C_L_endurance_UAS/a)+alpha_L_0_UAS;
v_endurance_aircraft = sqrt(2*GTOW/(C_L_endurance_aircraft*rho*S_wing));
v_endurance_UAS = sqrt(2*GTOW/(C_L_endurance_UAS*rho*S_wing));
fprintf("Velocity that corresponds to max endurance for the whole aircraft is: %f\n",...
    v_endurance_aircraft);
fprintf("Velocity that corresponds to max endurance for the UAS is: %f\n",...
    v_endurance_UAS);
fprintf("Angle of attack for max endurance of the whole aircraft is: %f\n",...
    alpha_endurance_max_aircraft);
fprintf("Angle of attack for max endurance of the UAS is: %f\n",...
    alpha_endurance_max_UAS);
range = 8*L_D_aircraft_max;
fprintf("Attainable range for the aircraft is: %f\n",range);
%% Calculating h_np
h_np = .25 + .6*(.0625/.0870);