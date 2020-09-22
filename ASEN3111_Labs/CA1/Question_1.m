function [] = Question_1(const)
%Q1_TRAPEZOIDAL Performs all calculations and outputs for question 1
%
% Author: Johnathan Tucker
% Collaborators: N/A
% This function takes in the struct of constants and outputs all of the
% required command line statements and plots for the first question
%
% Last Revised: 2/2/2020
%% First obtain the exact solution using MATLABs symbolic solver
% Notice that since there is no AoA N=L and A=D
syms P(theta)
% Define pressure as a function of theta
P(theta) = const.q_inf*(-4*(sin(theta))^2 -4*sin(theta)) + const.p_inf;
% Define L and D and Solve them symbolically
L_exact = double(-const.r*int(P*sin(theta),theta,0,2*pi));
D_exact = double(-const.r*int(P*cos(theta),theta,0,2*pi));
fprintf("%%%%%%%%%%%%%%%%BEGIN QUESTION ONE ANSWERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
fprintf("Exact Lift Solution is: %d [N/m]\n",L_exact)
fprintf("Exact Drag Solution is: %d [N/m]\n\n",D_exact)

%% Apply trapezoid rule to get lift and drag
% Pre-allocate memory
D_trap_vec = zeros(10,1);
L_trap_vec = zeros(10,1); %
% Loop through the trapezoid rule 9 times
for k = 2:10
    % Increment through the number of panels
    num_panels = linspace(0,2*pi,k);
    % Reset sum variables each loop
    sum_var_L = 0;
    sum_var_D = 0;
    for i = 1:k-1
        % Perform next trapezoidal iteration for lift
        sum_var_L_1 = 0.5*abs(num_panels(i+1) - num_panels(i))*...
            double((P(num_panels(i+1))*sin(num_panels(i+1)) +...
            P(num_panels(i))*sin(num_panels(i))));
        % Add it to the ongoing sum
        sum_var_L = sum_var_L + sum_var_L_1;
        % Update the lift solution
        L_trap = -const.r*sum_var_L;
        
        % Perform next trapezoidal iteration for drag
        sum_var_D_1 = 0.5*abs(num_panels(i+1) - num_panels(i))*...
            double((P(num_panels(i+1))*cos(num_panels(i+1)) +...
            P(num_panels(i))*cos(num_panels(i))));
        % Add it to the ongoing sum
        sum_var_D = sum_var_D + sum_var_D_1;
        % Update the drag solution
        D_trap = -const.r*sum_var_D;
    end
    % Save obtained L and D into their respective vectors
    L_trap_vec(k,:) = L_trap;
    D_trap_vec(k,:) = D_trap;
end

% Get the absolute error between the exact and estimated values
trap_L_error = abs(L_trap_vec(end) - L_exact);
trap_D_error = abs(D_trap_vec(end) - D_exact);

% These panel values were obtained from plot inspection
L_trap_panels = 4;
D_trap_panels = 3;

% Print out necessary information for the trapezoidal solutions
fprintf("Trapezoidal Lift solution is: %d [N/m]\n",L_trap_vec(end))
fprintf("Error between trapezoidal lift solution and exact lift solution is %d\n",trap_L_error)
fprintf("It took %d panels to get trapezoidal Lift within 0.001 N of the exact Lift\n\n",L_trap_panels)

fprintf("Trapezoidal Drag solution is: %d [N/m]\n",D_trap_vec(end))
fprintf("Error between trapezoidal drag solution and exact drag solution is %d\n",trap_D_error)
fprintf("It took %d panels to get trapezoidal Drag within 0.001 N of the exact Drag\n\n",D_trap_panels)

%% Apply Simpsons rule to get lift and drag
% Loop through the simpsons rule process 10 times 
for i = 1:10
    % Recalculate and reassign necessary variables each loop
    h = 2*pi/(2*i);
    sum_var = 0;
    % Loop through the simpsons process for drag
    for k = 1:i
        % Calculate the next increment in the summation
        sum_var_1 = double(P((2*k-2)*h)*cos((2*k-2)*h) +...
        4*P((2*k-1)*h)*cos((2*k-1)*h) + P((2*k)*h)*cos((2*k)*h));
        % Add it to all previous increments
        sum_var = sum_var + sum_var_1;
        % Update drag solution
        D_simps = -(h/3)*const.r*sum_var;

    end
    % Re-zero summation variable
    sum_var = 0;
    % Loop through the simpsons process for lift
    for k = 1:i
        % Calculate the next increment in the summation
        sum_var_1 = double(P((2*k-2)*h)*sin((2*k-2)*h) +...
        4*P((2*k-1)*h)*sin((2*k-1)*h) + P((2*k)*h)*sin((2*k)*h));
        % Add it to all previous increments
        sum_var = sum_var + sum_var_1;
        % Update drag solution
        L_simps = -(h/3)*const.r*sum_var;

    end
    % Save each increment to a vector for plotting
    L_simps_vec(i,:) = L_simps;
    D_simps_vec(i,:) = D_simps;
end
% Calculate absolute error between simpsons solution and exact solution
simp_L_error = abs(L_simps_vec(end) - L_exact);
simp_D_error = abs(D_simps_vec(end) - D_exact);

% These panel values were obtained from plot inspection
L_simp_panels = 3;
D_simp_panels = 2;
%% Print out necessary information for the simpsons solutions
fprintf("Simpsons Lift solution is: %d [N/m]\n",L_simps_vec(end))
fprintf("Error between simpsons lift solution and exact lift solution is: %d\n",simp_L_error)
fprintf("It took %d panels to get simpson Lift within 0.001 N of the exact Lift\n\n",L_simp_panels)

fprintf("Simpsons Drag solution is: %d [N/m]\n",D_simps_vec(end))
fprintf("Error between simpsons drag solution and exact drag solution is: %d\n",simp_D_error)
fprintf("It took %d panels to get simpson Drag within 0.001 N of the exact Drag\n\n",D_simp_panels)

% Create Plots for question 1
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plot(1:1:10, L_trap_vec)
hold on
plot(1:1:10, L_simps_vec)
xlabel("$Number\:of\:Panels$",'Interpreter','latex','FontSize',18)
ylabel("$Lift\:per\:Unit\:Span[N/m]$",'Interpreter','latex','FontSize',18)
title("$Simpsons\:and\:Trapezoidal\:Lift\:per\:Unit\:Span\:vs.\:Number\:of\:Panels$",...
    'Interpreter','latex','FontSize',18)
legend("$Trapezoidal\:Method$","$Simpsons\:Method$",'Interpreter','latex')

figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plot(1:1:10, D_trap_vec)
hold on
plot(1:1:10, D_simps_vec)
xlabel("$Number\:of\:Panels$",'Interpreter','latex','FontSize',18)
ylabel("$Drag\:per\:Unit\:Span[N/m]$",'Interpreter','latex','FontSize',18)
title("$Simpsons\:and\:Trapezoidal\:Drag\:per\:Unit\:Span\:vs.\:Number\:of\:Panels$",...
    'Interpreter','latex','FontSize',18)
legend("$Trapezoidal\:Method$","$Simpsons\:Method$",'Interpreter','latex')

end

