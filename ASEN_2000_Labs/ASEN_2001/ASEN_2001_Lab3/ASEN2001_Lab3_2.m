%% Housekeeping
clear; clc; close all;
%% Solving for p0 symbolically
% Declare symbolic variables
syms x p0
% Create constants
L = 36/39.37;
t_balsa = 0.03125/39.37;
t_foam = 0.75/39.37;
t_total = 2*t_balsa + t_foam;
h_f = 0.75 / 39.37;
h_b = (1/32) / 39.37;
c = (h_f + 2* h_b) / 2;
[I_b, I_f] = beammoment(4/39.37);
foam_modulus = 0.035483*(10^9);
balsa_modulus = 3.2953*(10^9);
%Get support reactions
R = int((4/39.37)*p0*sqrt(1-(2*x/L)^2),x,-L/2,L/2)/2;
%Get shear force
V = R - int((4/39.37)*p0*sqrt(1-(2*x/L)^2),x,-L/2,x);
%Calculate moment
M = int(V,x,-L/2,x);
sigma = (1.6174*10^7)/1.5;
right_side = (sigma*(I_b + (foam_modulus/balsa_modulus)*I_f))/c;
%Solve for p0 and save to variable another
M_something = subs(M,x,0);
another = double(solve(M_something == right_side,p0));
%Sub the 'another' variable into p0 for moment
M_real = subs(M,p0,another);
x_grid = linspace(-L/2,L/2,100);
M_plot = double(subs(M_real,x,x_grid));

%Moment plot
figure(1)
p(1) = plot(x_grid,M_plot,'DisplayName','Moment');
hold on

shear = (7.3894*10^4)/1.2;
A_f = (4/39.37)*t_foam;
right_side_2 = shear*(3/2)/A_f;
% Sub the 'another' variable into p0 for shear 
V_p0 = subs(V,p0,another);
V_plot = abs(double(subs(V_p0,x,x_grid)));
%Shear plot
p(2) = plot(x_grid,V_plot,'DisplayName','Shear');
title('$Shear\:and\:Moment\:Diagrams$','Interpreter','latex')
xlabel('$x\:(m)$','Interpreter','latex')
ylabel('$Shear\:(N)\:and\:Moment\:(Nm)$','Interpreter','latex')
legend([p(1),p(2)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving for the width symbolically
syms width
I_f = (1/12)*width*(h_f^3);
I_b = 2 * ((1/12) * width * (h_b)^3 + (width * h_b) * (0.5 * (h_f + h_b))^2);
eqn = (c*M_real)/(I_b + (foam_modulus/balsa_modulus)*I_f) == sigma;
%This should be the bending stress width function
w_moment = solve(eqn, width);
%Plot it by substituting x values in
w_plot = double(subs(w_moment,x,x_grid))/2;

figure(2)
plot(x_grid,w_plot)
hold on

A_f = width*t_foam;
%This should be the shear stress width function
w_shear = solve(((3/2)*V_p0/A_f)== shear,width);
w_shear_plot = abs(double(subs(w_shear,x,x_grid)))/2;

plot(x_grid,w_shear_plot)
title('$Shear\:and\:Moment\:Width\:Functions$','Interpreter','latex')
xlabel('$x\:(m)$','Interpreter','latex')
ylabel('$Width\:(m)$','Interpreter','latex')
hold on

% Find maximum force the wing will support
syms x
q = (4*0.0254) * another * sqrt(1-(2*x/L)^2);
max_load = subs(int(q,x,[-L/2 L/2]));
fprintf('The maximum force the wing will support is %0.3f N.\n', double(max_load))
%% Finding the centroid locations and wiffle brace locations
% [x_locations,centroids,distances,tier_2] = calc_x_locations(another);
% hold on
% vline(centroids(1))
% hold on
% vline(centroids(2))
% hold on
% vline(centroids(3))
% hold on
% vline(centroids(4))
% hold on
% vline(centroids(5))
% hold on
% vline(centroids(6))
% hold on
% vline(centroids(7))
% hold on
% vline(centroids(8))
% 
% vector = [-.4572, -.3429, -.2286,-.1143,0,.1143,.2286,.3429,.4572];
% hold on
% vline(vector(1),'b')
% hold on
% 
% vline(vector(2),'b')
% hold on
% vline(vector(3),'b')
% hold on
% vline(vector(4),'b')
% hold on
% vline(vector(5),'b')
% hold on
% vline(vector(6),'b')
% hold on
% vline(vector(7),'b')
% hold on
% vline(vector(8),'b')
% hold on
% vline(vector(9),'b')
% % 
% % hold on
% % vline(x_locations(1),'g')
% % hold on
% % vline(x_locations(2),'g')
% % hold on
% % vline(x_locations(3),'g')
