function [] = Question_2(const)
%QUESTION_2 Performs all calculations and outputs for question 2
%
% Author: Johnathan Tucker
% Collaborators: N/A
% This function takes in the struct of constants and outputs all of the
% required command line statements and plots for the second question
%
% Last Revised: 2/2/2020
%% Begin the solution process for question 2
% Load the Cp data
load('Cp.mat')

% Solve for lift using the trapezoid method. This will be used for the
% relative error when determining the number of panels for 1,5, and .1
% percent error.

%NOTE: The code commented below was run through the trapezoidal rule with
%150,000 panels to get an approximate lift value that will be used as the
%"exact" value. Uncomment the code below to confirm the L_exact_2 and
%D_exact_2 values below it.
% num_panels = 150000;
% [L_exact_2, D_exact_2] = airfoil_LD(const,num_panels,Cp_upper,Cp_lower);
L_exact_2 = 1.188109731345318e+03; %N/m
D_exact_2 = 1.46345878973628; %N/m
fprintf("%%%%%%%%%%%%%%%%BEGIN QUESTION TWO ANSWERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
fprintf("The Lift per unit span over the airfoil was caluclated to be %d [N/m] using the trapezoidal rule.\n",L_exact_2)
fprintf("The Drag per unit span over the airfoil was caluclated to be %d [N/m] using the trapezoidal rule.\n\n",D_exact_2)

% Create constants necessary to begin while loops
L_error = 0;
num_panels = 1;
% Loop until the relative error is equal to or less than 5 percent
while abs(L_exact_2 - L_error)/L_exact_2 > 0.05
    % Get the updated Lift
    [L_error,~] = airfoil_LD(const,num_panels,Cp_upper,Cp_lower);
    % Update the 5 percent number of panels solution
    num_panels_5_percent = num_panels;
    % Increment number of panels
    num_panels = num_panels + 1;
end
% Loop until the relative error is equal to or less than 1 percent
while abs(L_exact_2 - L_error)/L_exact_2 > 0.01
    % Get the updated Lift
    [L_error,~] = airfoil_LD(const,num_panels,Cp_upper,Cp_lower);
    % Update the 1 percent number of panels solution
    num_panels_1_percent = num_panels;
    % Increment number of panels
    num_panels = num_panels + 1;
end
% Seed the number of panels to 600 to reduce runtime. This value is less
% than the actual value but determined by letting it run completely without
% the seed a few times.
num_panels = 960;
% Loop until the relative error is equal to or less than .1 percent
while abs(L_exact_2 - L_error)/L_exact_2 > 0.001
    % Get the updated Lift
    [L_error,~] = airfoil_LD(const,num_panels,Cp_upper,Cp_lower);
    % Update the .1 percent number of panels solution
    num_panels_tenth_percent = num_panels;
    % Increment number of panels
    num_panels = num_panels + 1;
end
%% Begin Required Outputs Section
% Print all required statements
fprintf("The number of panels required for a 5%% relative error lift solution is: %d\n",num_panels_5_percent);
fprintf("The number of panels required for a 1%% relative error lift solution is: %d\n",num_panels_1_percent);
fprintf("The number of panels required for a 0.1%% relative error lift solution is: %d\n",num_panels_tenth_percent);

% The chunk of code below will compute the lift for a scaling number of
% panels up to 1500. The lift value at each panel was put into a vector and
% that vector was then saved to be used later.
% for i = 1:1500
%     [L,D] = airfoil_LD(const,i,Cp_upper,Cp_lower);
%     L_vec(i,:) = L;
% end
% save('LiftVector','L_vec')


% Convert number of panels to number of integration points
N = 1:1:1500;
n = N + 1;
% Use the above saved lift vector to plot the lift vs number of integration
% steps n
load('LiftVector.mat')
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plot(n,L_vec)
hold on
xline(num_panels_5_percent+1,'r');
hold on
xline(num_panels_1_percent+1, 'g');
hold on
xline(num_panels_tenth_percent+1, 'b');
xlabel("$Number\:of\:Integration\:Points\:n$",'Interpreter','latex','FontSize',18)
ylabel("$Lift\:per\:Unit\:Span[N/m]$",'Interpreter','latex','FontSize',18)
title("$Lift\:per\:Unit\:Span\:vs.\:Number\:of\:Integration\:Points$",...
    'Interpreter','latex','FontSize',18)
legend("Lift","$5\%\:Relative\:Error$","$1\%\:Relative\:Error$","$0.1\%\:Relative\:Error$",'Interpreter','latex')
end

