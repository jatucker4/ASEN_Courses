function Question_2()
%Question_2 Performs all calculations and outputs for question 2 in CA3
%
% Author: Johnathan Tucker
%
% Collaborators: N/A
%
% This function has no inputs or direct outputs. However, it does display
% the number of panels required to achieve a CL of 0.01 to the terminal. In
% addition to this terminal display it also displays three plot. The first
% is the convergence of CL as the number of panels increases. The second is
% the CL of the NACA 0012 at different AOA's. The final is the CP of the
% NACA 0012 at different AOA's.
% Last Revised: 3/26/2020
%% Begin code block for problem 2
% To determine a measurement of panels need to aquire accuracy I'll use the
% CP values from the vortex panel method to calculate CL and compare these
% with what should be the exact value for the NACA 0012 at zero degrees
% AOA, CL = 0;

% Exact CL value const
CL_0012_exact = 0;

% Now iterate through different panel values noting when 1 percent error is
% achieved
error_cl = 100;
% Get the panel increment close so it doesn't run forever
panel_increment = 927;
% Iterate on the process until CL is 0.01 or lower
while error_cl >= 0.01
    % Iterate on calculating the x and y values and then CL and CP
    [naca_0012_x_iter,naca_0012_y_iter] = NACA_Airfoil(0/100,0/10,12/100,1,panel_increment);
    [Cl_0012_iter,Cp_0012_iter] = Vortex_Panel(naca_0012_x_iter,naca_0012_y_iter,1,0,0,0);
    % Increment the panel number
    panel_increment = panel_increment + 1;
    % Calculate the error between the exact and calculated CL values
    error_cl = abs(CL_0012_exact - Cl_0012_iter);
end
fprintf("The number of panels required to obtain 0.01 absolute error...\nwhen compared to the exact value is %f\n\n",panel_increment)

%% Create a plot for CL versus the number of panels

% error_cl = 100;
% % Get the panel increment close so it doesn't run forever
% panel_increment = 4;
% % Iterate on the process until CL is 0.01 or lower
% while panel_increment ~= 929
%     % Iterate on calculating the x and y values and then CL and CP
%     [naca_0012_x_iter,naca_0012_y_iter] = NACA_Airfoil(0/100,0/10,12/100,1,panel_increment);
%     [Cl_0012_iter(panel_increment - 3),Cp_0012_iter] = Vortex_Panel(naca_0012_x_iter,naca_0012_y_iter,1,0,0,0);
%     % Increment the panel number
%     panel_increment = panel_increment + 1;
%     % Calculate the error between the exact and calculated CL values
%     error_cl = abs(CL_0012_exact - Cl_0012_iter);
% end

% Load in the full CL vector. This vector was created by running the code
% above that is commented out.
CL_iter_vec = load("CL_iteration_mat.mat");
% Save the CL vector to a const.
CL_iter_vec = CL_iter_vec.Cl_0012_iter;
% Create the Figure
figure
plot(4:929,[CL_iter_vec,CL_iter_vec(end)])
hold on
plot(1:928,ones(1,928).*0.01)
title('CL vs Number of Panels for NACA 0012',...
    'Interpreter','latex','FontSize',16)
xlabel('$Panel\:Number$','Interpreter','latex','FontSize',12)
ylabel('$Coefficient\:of\:Lift$','Interpreter','latex','FontSize',12)
legend('$CL$','$0.01\:Desired\:CL$','Interpreter','latex')

%% Now using the panel increment to reach the desired error I'll create
% plots of CP at the desired angles of attack 
alpha_vec = [-5,0,5,10];
figure
% Iterate through the desired AOA's and get the corresponding CL's
for i = 1:length(alpha_vec)
    [Cl_0012_iter_2(i),Cp_0012_iter] = Vortex_Panel(naca_0012_x_iter,naca_0012_y_iter,1,alpha_vec(i),2,i);
end
%% Finally create a plot of the CL values at these specified AOA's
figure
% Iterate through the CL values and plot them
for i = 1:length(Cl_0012_iter_2)
    scatter(alpha_vec(i),Cl_0012_iter_2(i))
    hold on
end
title('$Coefficeint\:of\:Lift\:vs\:\alpha\:for\:NACA\:0012$',...
    'Interpreter','latex','FontSize',18)
xlabel('$\alpha$','Interpreter','latex','FontSize',12)
ylabel('$Coefficient\:of\:Lift$','Interpreter','latex','FontSize',12)
legend('CL at $\alpha$=-5','CL at $\alpha$=0','CL at $\alpha$=5',...
    'CL at $\alpha$=10','Interpreter','latex')
end

