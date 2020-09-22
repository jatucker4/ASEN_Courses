%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When called this function conducts sensitivity analysis on method 3 and
% outputs a plot of the coefficient of restitutions sensitivity to varying
% amounts of uncertainty in its inputs assuming a standard bounce time
% vector.
%
% Created by: Johnathan Tucker
% 
% Inputs: figure_number: Figure number for organized plots
%
% Outputs: A plot of the sensitivity of the coefficient of restitution
%          sens_vec_t_total: A vector of the sensitivity of e to the time of
%                            the previous bounce in seconds
%          sens_vec_h_initial: A vector of the sensitivity of e to the current
%                              bounce time in inches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sens_vec_t_total,sens_vec_h0] = sens_method_3(figure_number)
% Create a vector of the total times in seconds
total_time_vec = [7.6900,8.8000,8.0500,8.7000,8.3900,7.2900,8.9700,7.1500,...
    8.4300,8.6700];
% Create a vector of the uncertainty in the times in seconds
time_unc_vec = linspace(0,1,10);
% Create the initial height value in inches
initial_height = 28;
% Create the height uncertainty in inches
unc_vec = linspace(0,1,10);
% Create a gravity variable in (in/s^2)
g = 386.09;
for i = 1:length(total_time_vec)
    % Calculate the uncertainty of e by varying only the uncertainty in the
    % total time
    sens_vec_t_total(i) = sqrt( (((total_time_vec(i) + sqrt((2 * ...
        initial_height)/g)) - ((total_time_vec(i) - sqrt((2 *...
        initial_height)/g))))/(((total_time_vec(i) + sqrt((2 *...
        initial_height)/g)))^2) * time_unc_vec(i))^2 + ...
        (((total_time_vec(i) + sqrt((2 *initial_height)/g))*(1/g)*...
        ((2*initial_height/g)^(-1/2)) + (total_time_vec(i) - ...
        sqrt((2 *initial_height)/g))*(1/g)*...
        ((2*initial_height)/g)^(-1/2))/(((total_time_vec(i) + sqrt((2 *...
        initial_height)/g)))^2) * 0.3)^2);
    % Calculate the uncertainty of e by varying only the uncertainty in the
    % initial height
    sens_vec_h0(i) = sqrt( (((total_time_vec(i) + sqrt((2 * ...
        initial_height)/g)) - ((total_time_vec(i) - sqrt((2 *...
        initial_height)/g))))/(((total_time_vec(i) + sqrt((2 *...
        initial_height)/g)))^2) * 0.05)^2 + ...
        (((total_time_vec(i) + sqrt((2 *initial_height)/g))*(1/g)*...
        ((2*initial_height/g)^(-1/2)) + (total_time_vec(i) - ...
        sqrt((2 *initial_height)/g))*(1/g)*...
        ((2*initial_height/g)^(-1/2)))/(((total_time_vec(i) + sqrt((2 *...
        initial_height)/g)))^2) * unc_vec(i))^2);
end
% Create the figure
figure(figure_number)
% Plot the sensitivity of e to the initial height
plot(unc_vec,sens_vec_h0,'LineWidth',2)
title('$Sensitivity\:of\:Coefficient\:of\:Restitution\:to\:Initial\:Height\:in\:Method\:3$',...
    'Interpreter','latex',...
    'FontSize',26)
xlabel('$Initial\:Height\:Uncertainty(in)$','Interpreter','latex','FontSize',26)
ylabel('$Coefficient\:of\:Restitution\:Sensitivity$','Interpreter','latex','FontSize',26)
legend('Initial Bounce Height')

% Create the figure
figure(figure_number+1)
% Plot the sensitivity of e to the total time
plot(unc_vec,sens_vec_t_total,'LineWidth',2)
title('$Sensitivity\:of\:Coefficient\:of\:Restitution\:to\:Total\:Time\:in\:Method\:3$',...
    'Interpreter','latex',...
    'FontSize',26)
xlabel('$Total\:Time\:Uncertainty(s)$','Interpreter','latex','FontSize',26)
ylabel('$Coefficient\:of\:Restitution\:Sensitivity$','Interpreter','latex','FontSize',26)
legend('Total Time')

end