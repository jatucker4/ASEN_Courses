function [arc_length, pos_final, time, z_t_vec, ...
    g_vec] = Lab1_Hill(coaster_plot, pos_initial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by: Johnathan Tucker in collaboration with Elliot Tung
% Purpose: Create top hat piece for the coaster track and output
%          necessary information for other pieces and calculations
% Inputs: 
%         pos_initial: The x,y,z coordinate of the beginning of the banked
%                    turn in meters
%         coaster_plot: The figure name to plot the banked portion
%
% Outputs: pos_final: Final x,y,z coordinates at the end point of the
%                     banked turn in meters
%          arc_length: Length of the banked turn meters
%          time: Time it took to complete the banked turn in seconds
%          z_t_vec: The change in height over time of the banked turn in
%                   meters
%          g_vec: Vector of G's felt during the top hat portion. Unitless
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parametric Function of Coaster
x_initial = pos_initial(1);
y_initial = pos_initial(2);
z_initial = pos_initial(3);
tot_length = 0;
% Crete first segment parametric equation
syms t k
x1 = -20*cos(t) + x_initial;
y1 = -0*t + y_initial;
z1 = 10+20*sin(t) + z_initial;
s_t = matlabFunction(sqrt((diff(x1)^2) + (diff(z1)^2)));
arc_length = integral(s_t ,3*pi/2,2*pi);
tot_length = tot_length + arc_length;
% Calculate G's for first section
g_1 = 2*(125 - k)/20 + 1;
tspan = linspace(3*pi/2,2*pi);
z_t_vec_1 = double(subs(z1,t,tspan));
g_vec_1 = double(subs(g_1, k, z_t_vec_1));
% Plot first section
figure(coaster_plot)
fplot3(x1,y1,z1,[3*pi/2,2*pi], 'LineWidth',2)
hold on
scatter3(double(subs(x1,t,2*pi)), double(subs(y1,t,2*pi)), double(subs(z1,t,2*pi)), 'filled')
hold on
axis equal
% Create second segment parametric equation
syms t
x2 = 10 - 0*t + 82.29;
y2 = -0*t + 72.99;
z2 = t - 1.98 + 2.01;
tot_length = tot_length + (114.9 - 85);
% Calculate G's for second section
tspan = linspace(95,114.9);
z_t_vec_2 = double(subs(z2,t,tspan));
g_vec_2 = zeros(1,length(z_t_vec_2));
% Plot second section
fplot3(x2,y2,z2,[95,114.9], 'LineWidth',2)
hold on
scatter3(double(subs(x2,t,114.9)), double(subs(y2,t,114.9)), double(subs(z2,t,114.9)), 'filled')
hold on
% Create third segment parametric equation
x3 = 10 - 0*t + 82.29;
y3 = -10*cos(t)-10 + 72.99;
z3 = 114.9 + 10*sin(t);
s_t = matlabFunction(sqrt((diff(y3)^2) + (diff(z3)^2)));
arc_length = integral(s_t ,0,pi);
tot_length = tot_length + arc_length;
% Calculate g's for third segment
tspan = linspace(0,pi);
z_t_vec_3 = double(subs(z3,t,tspan));
g_vec_3 = double(subs(g_1,k,z_t_vec_3));
% Plot third section
fplot3(x3,y3,z3,[0,pi], 'LineWidth',2)
hold on
scatter3(double(subs(x3,t,pi)), double(subs(y3,t,pi)), double(subs(z3,t,pi)), 'filled')
hold on
% Create fourth segment parametric equations
x4 = 10 - 0*t + 86.93 - 4.64;
y4 = -0*t - 20 + .99 +19.01 + 52.99;
z4 = t;
tot_length = tot_length + (114.9 - 40);
% Calculate g's of fourth section
tspan = linspace(114.9,40);
z_t_vec_4 = double(subs(z4,t,tspan));
g_vec_4 = zeros(1,length(z_t_vec_4));
% Plot fourth section 
fplot3(x4,y4,z4,[40,114.9], 'LineWidth',2)
hold on
scatter3(double(subs(x4,t,114.9)), double(subs(y4,t,114.9)), ...
    double(subs(z4,t,114.9)), 'filled')
hold on
% Create fifth section parametric equations
x5 = -40*cos(t) -30 + 86.93 + 75.36;
y5 = -0*t - 20 + .99 +19.01 + 52.99;
z5 = 40+40*sin(t);
s_t = matlabFunction(sqrt((diff(x5)^2) + (diff(z5)^2)));
arc_length = integral(s_t ,3*pi/2,2*pi);
tot_length = tot_length + arc_length;
% Calculate the g's for fifth section
tspan = linspace(2*pi,3*pi/2);
z_t_vec_5 = double(subs(z5,t,tspan));
g_2 = 2*(125 - k)/40 - 1;
g_vec_5 = double(subs(g_2,k,z_t_vec_5));
% Plot fifth segment
fplot3(x5,y5,z5,[3*pi/2,2*pi], 'LineWidth',2)
hold on
scatter3(double(subs(x5,t,2*pi)), double(subs(y5,t,2*pi)), ...
    double(subs(z5,t,2*pi)), 'filled')
arc_length = tot_length;
pos_final = [double(subs(x5,t,3*pi/2)),double(subs(y5,t,3*pi/2)),...
    double(subs(z5,t,3*pi/2))];
hold on
scatter3(double(subs(x5,t,3*pi/2)),double(subs(y5,t,3*pi/2)),...
    double(subs(z5,t,3*pi/2)), 'filled')
time = 10;
% Create output vectors
z_t_vec = [z_t_vec_1, z_t_vec_2, z_t_vec_3, z_t_vec_4, z_t_vec_5];
g_vec = [g_vec_1, g_vec_2, g_vec_3, g_vec_4, g_vec_5];
end