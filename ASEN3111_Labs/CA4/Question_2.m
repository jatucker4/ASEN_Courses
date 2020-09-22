function Question_2()
%Question_2 Performs all calculations and outputs for question 2 in CA4
%
% Author: Johnathan Tucker
%
% Collaborators: N/A
%
% This function has no inputs or direct outputs. However, it does display
% the number of odd numbers required to achieve 5, 1, and 0.1 percent relative
% error between an "exact" Lift and induced drag and the calculatedd Lift
% and induced drag.
%
% Last Revised: 4/8/2020
%% Create constants
b = 100; % [ft]
c_r = 15; % [ft]
c_t = 5; % [ft]
geo_r = 5; % [deg]
geo_t = 0; % [deg]
V_inf = 150*5280/3600; % [ft/s]
rho_SL = 0.0023769; % [slugs/ft^3]
S = (c_r + c_t)*(b/2); % [ft^2]
%% Get the root and tip aerodynamic twist
% First get the x and y values for each airfoil
[naca_0012_x,naca_0012_y] = NACA_Airfoil(0/100,0/10,12/100,1,150);
[naca_2412_x,naca_2412_y] = NACA_Airfoil(2/100,4/10,12/100,1,150);

% Then get the coefficient of lift for the airfoils for each angle of
% attack
aoa_vec = linspace(-5,10);
for i = 1:length(aoa_vec)
    [Cl_0012(i),~] = Vortex_Panel(naca_0012_x,naca_0012_y,V_inf,aoa_vec(i),0,0);
    [Cl_2412(i),~] = Vortex_Panel(naca_2412_x,naca_2412_y,V_inf,aoa_vec(i),0,0);
end

% Now get aero_r and a0_r
fit_2412 = polyfit(aoa_vec, Cl_2412,1);
aero_r = fit_2412(2)/fit_2412(1);
a0_r = fit_2412(1)*180/pi;

% Now get aero_r and a0_r
fit_0012 = polyfit(aoa_vec, Cl_0012,1);
aero_t = -fit_0012(2)/fit_0012(1);
a0_t = fit_0012(1)*180/pi;

%% Begin error calculations
% Now get an "exact" Lift and induced drag using a high number of panels
N = 1000;
[~,exact_cl,exact_cdi] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
exact_L = 0.5*rho_SL*(V_inf^2)*S*exact_cl;
exact_Di = 0.5*rho_SL*(V_inf^2)*S*exact_cdi;

exact_L_SI = 0.5*rho_SL*(V_inf^2)*S*exact_cl*4.44822;
exact_Di_SI = 0.5*rho_SL*(V_inf^2)*S*exact_cdi*4.44822;
% Print out the calculated lift and drag in Newtons
fprintf("For 1000 odd terms at sea level and a velocity of 150 mph:\n");
fprintf("Lift = %f [N]\n",exact_L_SI)
fprintf("Induced Drag = %f [N]\n\n",exact_Di_SI)
% Using the exact values calculate the number of panels it takes to achieve
% 5, 1, and 0.1 percent relative error
N = 1;
error_L = 100;
error_Di = 100;

% Create a loop to get the 0.1 percent error
while true
    [~,cl,cdi] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
    L = 0.5*rho_SL*(V_inf^2)*S*cl;
    Di = 0.5*rho_SL*(V_inf^2)*S*cdi;
    
    error_L(N) = (abs(L - exact_L)/exact_L) * 100;
    error_Di(N) = (abs(Di - exact_Di)/exact_Di) * 100;
    
    if error_L(N) <=0.1 && error_Di(N) <= 0.1
        break
    end
    N = N + 1;
end

% Now find the number of panels where the error is first below 5 percent
% for both Lift and Drag
num_five_perc_L = find(lt(error_L,5),1,'first');
num_five_perc_Di = find(lt(error_Di,5),1,'first');
% Print the number of panels need for five percent error
fprintf("The number of odd terms required for a 5%% relative error in \nthe lift and induced drag solutions is: %d\n",num_five_perc_L);
fprintf("The number of odd terms required for a 5%% relative error in \nthe lift and induced drag solutions is: %d\n\n",num_five_perc_Di);

% Now find the number of panels where the error is first below 1 percent
% for both Lift and Drag
num_one_perc_L = find(lt(error_L,1),1,'first');
num_one_perc_Di = find(lt(error_Di,1),1,'first');
% Print the number of panels need for five percent error
fprintf("The number of odd terms required for a 1%% relative error in \nthe lift solution is: %d\n",num_one_perc_L);
fprintf("The number of odd terms required for a 1%% relative error in \nthe induced drag solution is: %d\n\n",num_one_perc_Di);


% Finally display the number of panels needed to achieve a relative error
% less than 0.1 percent in both lift and drag
num_point_1_perc_L = find(lt(error_L,.1),1,'first');
fprintf("The number of odd terms required for a 0.1%% relative error in \nthe lift solutions is: %d\n",num_point_1_perc_L);
fprintf("The number of odd terms required for a 0.1%% relative error in \nthe induced drag solutions is: %d\n\n",N);
end

