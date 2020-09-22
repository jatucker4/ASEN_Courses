function Question_3()
%Question_3 Performs all calculations and outputs for question 3 in CA3
%
% Author: Johnathan Tucker
%
% Collaborators: Michael Martinson
%
% This function has no inputs or direct outputs. However, it does display
% the error between thick and thin airfoil calculated zero lift angle of
% attack and lift curve slope to the terminal. In addition to terminal
% displays it also displays a plot of CL vs alpha for different NACA
% airfoils.
%
% Last Revised: 3/27/2020
%% Begin code block for problem 3
% First get the x and y coords for each airfoil
[naca_0012_x,naca_0012_y] = NACA_Airfoil(0/100,0/10,12/100,1,150);
[naca_2412_x,naca_2412_y] = NACA_Airfoil(2/100,4/10,12/100,1,150);
[naca_4412_x,naca_4412_y] = NACA_Airfoil(4/100,4/10,12/100,1,150);
[naca_2424_x,naca_2424_y] = NACA_Airfoil(2/100,4/10,24/100,1,150);

% Now get the Cl for each airfoil iterating alpha from -5 to 10 degrees
alpha_vec = linspace(-5,10,20);
for i = 1:length(alpha_vec)
    [Cl_0012(i),~] = Vortex_Panel(naca_0012_x,naca_0012_y,1,alpha_vec(i),0,0);
    [Cl_2412(i),~] = Vortex_Panel(naca_2412_x,naca_2412_y,1,alpha_vec(i),0,0);
    [Cl_4412(i),~] = Vortex_Panel(naca_4412_x,naca_4412_y,1,alpha_vec(i),0,0);
    [Cl_2424(i),~] = Vortex_Panel(naca_2424_x,naca_2424_y,1,alpha_vec(i),0,0);
end

% Plot the Cl for each airfoil on the same plot
figure
plot(alpha_vec,Cl_0012)
hold on
plot(alpha_vec,Cl_2412)
hold on
plot(alpha_vec,Cl_4412)
hold on
plot(alpha_vec,Cl_2424)
ylabel('$CL$','Interpreter','latex','FontSize',12)
xlabel('$\alpha$','Interpreter','latex','FontSize',12)
title('CL vs $\alpha$ of Various NACA Airfoils','Interpreter','latex','FontSize',18)
legend("$NACA\:0012$","$NACA\:2412$","$NACA\:4412$","$NACA\:2424$",...
    'Interpreter','latex')

% Compare each Cl with the thin airfoil theory Cl and the zero lift angle
% of attack with the thin airfoil theory version
% First Use polyfit to get the slope
slope_0012 = polyfit(alpha_vec.*pi/180,Cl_0012,1);
slope_2412 = polyfit(alpha_vec.*pi/180,Cl_2412,1);
slope_4412 = polyfit(alpha_vec.*pi/180,Cl_4412,1);
slope_2424 = polyfit(alpha_vec.*pi/180,Cl_2424,1);

% Calculate the percent difference between thick and thin airfoil theories
slope_error_0012 = 100 * abs(slope_0012(1) - (2*pi))/slope_0012(1);
slope_error_2412 = 100 * abs(slope_2412(1) - (2*pi))/slope_2412(1);
slope_error_4412 = 100 * abs(slope_4412(1) - (2*pi))/slope_4412(1);
slope_error_2424 = 100 * abs(slope_2424(1) - (2*pi))/slope_2424(1);

% Now get the zero lift angle of attack for the thick airfoil theory case
zeroL_aoa_0012 = -(180/pi)*slope_0012(2)/slope_0012(1);
zeroL_aoa_2412 = (180/pi)*slope_2412(2)/slope_2412(1);
zeroL_aoa_4412 = (180/pi)*slope_4412(2)/slope_4412(1);
zeroL_aoa_2424 = (180/pi)*slope_2424(2)/slope_2424(1);

% Get the zero lift angle of attack for the thin airfoil theory case
tat_zeroL_aoa_0012 = get_zlaoa_tat(0/100,0/10,1,150);
tat_zeroL_aoa_2412 = get_zlaoa_tat(2/100,4/10,1,150);
tat_zeroL_aoa_4412 = get_zlaoa_tat(4/100,4/10,1,150);
tat_zeroL_aoa_2424 = get_zlaoa_tat(2/100,4/10,1,150);

% Compute percent difference between the zero lift AOA values for each case
zeroL_error_0012 = abs(abs(zeroL_aoa_0012) - abs(tat_zeroL_aoa_0012));
zeroL_error_2412 = abs(abs(zeroL_aoa_2412) - abs(tat_zeroL_aoa_2412));
zeroL_error_4412 = abs(abs(zeroL_aoa_4412) - abs(tat_zeroL_aoa_4412));
zeroL_error_2424 = abs(abs(zeroL_aoa_2424) - abs(tat_zeroL_aoa_2424));

% Display all results for each airfoil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NACA 0012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NACA 0012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
fprintf("Results for the NACA 0012 Airfoil:\n")
fprintf("Vortex Panel Estimated Lift Curve Slope: %f [1/rad]\n",slope_0012(1))
fprintf("Thin Airfoil Theory Estimated Lift Curve Slope: %f [1/rad]\n",2*pi)
fprintf("Percent Error of Lift Curve Slope Estimates: %f %%\n\n",slope_error_0012)

fprintf("Vortex Panel Estimated Zero Lift AOA: %f [deg]\n",zeroL_aoa_0012)
fprintf("Thin Airfoil Theory Estimated Zero Lift AOA: %f [deg]\n",tat_zeroL_aoa_0012)
fprintf("Absolute Error of Zero Lift AOA Estimates: %f \n\n",zeroL_error_0012)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NACA 2412 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NACA 2412 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
fprintf("Results for the NACA 2412 Airfoil:\n")
fprintf("Vortex Panel Estimated Lift Curve Slope: %f [1/rad]\n",slope_2412(1))
fprintf("Thin Airfoil Theory Estimated Lift Curve Slope: %f [1/rad]\n",2*pi)
fprintf("Percent Error of Lift Curve Slope Estimates: %f %%\n\n",slope_error_2412)

fprintf("Vortex Panel Estimated Zero Lift AOA: %f [deg]\n",zeroL_aoa_2412)
fprintf("Thin Airfoil Theory Estimated Zero Lift AOA: %f [deg]\n",tat_zeroL_aoa_2412)
fprintf("Absolute Error of Zero Lift AOA Estimates: %f \n\n",zeroL_error_2412)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NACA 4412 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NACA 4412 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
fprintf("Results for the NACA 4412 Airfoil:\n")
fprintf("Vortex Panel Estimated Lift Curve Slope: %f [1/rad]\n",slope_4412(1))
fprintf("Thin Airfoil Theory Estimated Lift Curve Slope: %f [1/rad]\n",2*pi)
fprintf("Percent Error of Lift Curve Slope Estimates: %f %%\n\n",slope_error_4412)

fprintf("Vortex Panel Estimated Zero Lift AOA: %f [deg]\n",zeroL_aoa_4412)
fprintf("Thin Airfoil Theory Estimated Zero Lift AOA: %f [deg]\n",tat_zeroL_aoa_4412)
fprintf("Absolute Error of Zero Lift AOA Estimates: %f \n\n",zeroL_error_4412)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NACA 2424 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NACA 2424 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
fprintf("Results for the NACA 2424 Airfoil:\n")
fprintf("Vortex Panel Estimated Lift Curve Slope: %f [1/rad]\n",slope_2424(1))
fprintf("Thin Airfoil Theory Estimated Lift Curve Slope: %f [1/rad]\n",2*pi)
fprintf("Percent Error of Lift Curve Slope Estimates: %f %%\n\n",slope_error_2424)

fprintf("Vortex Panel Estimated Zero Lift AOA: %f [deg]\n",zeroL_aoa_2424)
fprintf("Thin Airfoil Theory Estimated Zero Lift AOA: %f [deg]\n",tat_zeroL_aoa_2424)
fprintf("Absolute Error of Zero Lift AOA Estimates: %f \n\n",zeroL_error_2424)
end

