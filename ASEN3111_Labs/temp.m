%% Main Script File for Comp Lab 4
% Author: Duncan McGough
% Created: 10/24/17

%% Housekeeping
clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2

% % This calls the NACA airfoil function for a NACA 0012 airfoil with 
% % 	c=1 and N=1000 points
% [x_exact y_exact] = NACA_Airfoil('0012', 1, 1000);
% 
% % Need to split the cp into an upper and lower curve
% [minx_exact minx_exact_index] = min(abs(x_exact));
% 
% % Split exact x locations into upper and lower values
% x_exact_l = x_exact(1:minx_exact_index);
% x_exact_u = x_exact(minx_exact_index+1:end);
% 
% % Find an accurate representation of Cp using N=1000. Velocity does not 
% % 	matter for comparison purposes. AoA = 0. 
% [cl_exact cp_exact] = Vortex_Panel(x_exact, y_exact, 10, 0, true);
% 
% % Split cp into upper and lower
% cp_exact_l = cp_exact(1:minx_exact_index);
% cp_exact_u = cp_exact(minx_exact_index+1:end);
% 
% % Integrate under cp curve for upper and lower
% cp_int_exact_l = trapz(cp_exact_l, x_exact_l);
% cp_int_exact_u = trapz(cp_exact_u, x_exact_u(1:end-1)); 
% % x is too long by one extra vector, so remove the redundant last point
% 
% 
% %% Start to iterate through number of panels
% for i=5:1000
% 	[x_test, y_test] = NACA_Airfoil('0012',1,i); % Calculate airfoil geometry
% 	[cl_test, cp_test] = Vortex_Panel(x_test, y_test, 10, 0, false); % create cp values
% 
% 
% 	%% Find where to split data into upper and lower
% 	[minx_test minx_test_index] = min(abs(x_test));
% 
% 	%% split x locations into upper and lower
% 	x_test_l = x_test(1:minx_test_index);
% 	x_test_u = x_test(minx_test_index+1:end);
% 
% 	%% split cp into upper and lower
% 	cp_test_l = cp_test(1:minx_test_index);
% 	cp_test_u = cp_test(minx_test_index+1:end);
% 
% 	%% Integrate
% 	cp_int_test_l = trapz(cp_test_l, x_test_l);
% 	cp_int_test_u = trapz(cp_test_u, x_test_u(1:end-1));
% 	% x is too long by one extra vector, so remove the redundant last point
% 
% 	%% Find error between test points and exact solution
% 	cp_error_l = abs( (cp_int_test_l - cp_int_exact_l)/(cp_int_exact_l) ); % lower error
% 	cp_error_u = abs( (cp_int_test_u - cp_int_exact_u)/(cp_int_exact_u) ); % upper error
% 
% 
% 	%% interpolate between the lower-resolution curve points for the higher resolution x values
% 	%cp_interp_u = interp1(x_test_u(1:end-1), cp_test_u, x_exact_u); 
% 	%cp_interp_l = interp1(x_test_l, cp_test_l, x_exact_l); 
% 	%cp_diff_u = abs(cp_interp_u-cp_exact_u); % find the difference between the exact value and the test value
% 	%cp_diff_l = abs(cp_interp_l-cp_exact_l); % find the difference between the exact value and the test value
% 	%cp_error_u = cp_diff_u./cp_exact_u; % find the associated error
% 	%cp_error_l = cp_diff_l./cp_exact_l; % find the associated error
% 
% 	%cp_error = mean([mean(cp_error_u), mean(cp_error_l)]);
% 
% 	% test to see what the error currently is
% 	if (cp_error_l <= 0.05) && (cp_error_u <= 0.05)
% 		n_iteration = i; % store the value of how many iterations it takes to get within a certain percentage
% 		%[cl_test, cp_test] = Vortex_Panel(x_test, y_test, 10, 0, true); % get plot for approximation
% 		break % if error is within 5% then break the for loop
% 	end
% 
% 
% 
% end
% 
% %% Plot several error plots
% [x_test1, y_test1] = NACA_Airfoil('0012',1,5); % Calculate airfoil geometry
% [cl_test1, cp_test1] = Vortex_Panel(x_test1, y_test1, 10, 0, true); % create cp values
% 
% [x_test2, y_test2] = NACA_Airfoil('0012',1,15); % Calculate airfoil geometry
% [cl_test2, cp_test2] = Vortex_Panel(x_test2, y_test2, 10, 0, true); % create cp values
% 
% [x_test3, y_test3] = NACA_Airfoil('0012',1,30); % Calculate airfoil geometry
% [cl_test3, cp_test3] = Vortex_Panel(x_test3, y_test3, 10, 0, true); % create cp values
% 
% %% Plot errorful plot versus accurate plot
% figure
% hold on
% grid on 
% h = gca;
% set(h,'YDir','reverse');
% set(h,'FontSize',18);
% plot(x_test(1:end-1),cp_test,'LineWidth',2)
% plot(x_exact(1:end-1), cp_exact, 'LineWidth',2)
% legend('Approximated Cp Curve', 'Exact Cp Curve')
% title('Cp vs Location')
% xlabel('x')
% ylabel('Cp')
% hold off
% 
% 
% %% Plot nominal panel data for different angles of attack
% [x_nom, y_nom] = NACA_Airfoil('0012',1,n_iteration);
% [cl_nom_n5, cp_nom_n5] = Vortex_Panel(x_nom, y_nom, 10, -5, false);
% [cl_nom_0, cp_nom_0] = Vortex_Panel(x_nom, y_nom, 10, 0, false);
% [cl_nom_5, cp_nom_5] = Vortex_Panel(x_nom, y_nom, 10, 5, false);
% [cl_nom_10, cp_nom_10] = Vortex_Panel(x_nom, y_nom, 10, 10, false);
% 
% %% Plot all plots 
% figure
% hold on
% grid on 
% h = gca;
% set(h,'YDir','reverse');
% set(h,'FontSize',18);
% plot(x_nom(1:end-1),cp_nom_n5,'LineWidth',2)
% plot(x_nom(1:end-1),cp_nom_0,'LineWidth',2)
% plot(x_nom(1:end-1),cp_nom_5,'LineWidth',2)
% plot(x_nom(1:end-1),cp_nom_10,'LineWidth',2)
% legend('AoA = -5*', 'AoA = 0*', 'AoA = 5*', 'AoA = 10*')
% title('Cp vs Location')
% xlabel('x')
% ylabel('Cp')
% hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 3

alpha_range = linspace(-5,10,30); % set up 50 alpha points to test

% Set up airfoil geometry for each NACA airfoil, using N=100 points
[x0012 y0012] = NACA_Airfoil('0012',1,150);
[x2212 y2212] = NACA_Airfoil('2212',1,150);
[x4412 y4412] = NACA_Airfoil('4412',1,150);
[x2430 y2430] = NACA_Airfoil('2430',1,150);

for i=1:length(alpha_range)
	aoa = (alpha_range(i));
	[cl0012(i) ~] = Vortex_Panel(x0012, y0012, 10, aoa, false);
	[cl2212(i) ~] = Vortex_Panel(x2212, y2212, 10, aoa, false);
	[cl4412(i) ~] = Vortex_Panel(x4412, y4412, 10, aoa, false);
	[cl2430(i) ~] = Vortex_Panel(x2430, y2430, 10, aoa, false);
end

% Find the lift slope 
slope_0012_rad = polyfit(alpha_range*pi/180, cl0012, 1);
slope_2212_rad = polyfit(alpha_range*pi/180, cl2212, 1);
slope_4412_rad = polyfit(alpha_range*pi/180, cl4412, 1);
slope_2430_rad = polyfit(alpha_range*pi/180, cl2430, 1);

% Find the lift slope for polyfit
slope_0012 = polyfit(alpha_range, cl0012, 1);
slope_2212 = polyfit(alpha_range, cl2212, 1);
slope_4412 = polyfit(alpha_range, cl4412, 1);
slope_2430 = polyfit(alpha_range, cl2430, 1);

% Find zero lift angle of attack
a0_0012 = -slope_0012(1,2)/slope_0012(1,1);
a0_2212 = -slope_2212(1,2)/slope_2212(1,1);
a0_4412 = -slope_4412(1,2)/slope_4412(1,1);
a0_2430 = -slope_2430(1,2)/slope_2430(1,1);
%% Plot all plots 

figure
hold on
grid on 
h = gca;
set(h,'FontSize',18);
plot(alpha_range, cl0012,'LineWidth',2)
plot(a0_0012,0,'o','LineWidth',2)
legend('Cl vs AoA', 'Zero-Lift AoA')
title('Cl vs AoA for a NACA 0012')
xlabel('AoA, degrees')
ylabel('Cl')
hold off

figure
hold on
grid on 
h = gca;
set(h,'FontSize',18);
plot(alpha_range, cl2212,'LineWidth',2)
plot(a0_2212,0,'o','LineWidth',2)
legend('Cl vs AoA', 'Zero-Lift AoA')
title('Cl vs AoA for a NACA 2212')
xlabel('AoA, degrees')
ylabel('Cl')
hold off

figure
hold on
grid on 
h = gca;
set(h,'FontSize',18);
plot(alpha_range, cl4412,'LineWidth',2)
plot(a0_4412,0,'o','LineWidth',2)
legend('Cl vs AoA', 'Zero-Lift AoA')
title('Cl vs AoA for a NACA 4412')
xlabel('AoA, degrees')
ylabel('Cl')
hold off

figure
hold on
grid on 
h = gca;
set(h,'FontSize',18);
plot(alpha_range, cl2430,'LineWidth',2)
plot(a0_2430,0,'o','LineWidth',2)
legend('Cl vs AoA', 'Zero-Lift AoA')
title('Cl vs AoA for a NACA 2430')
xlabel('AoA, degrees')
ylabel('Cl')
hold off
