%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When called this function conducts sensitivity analysis on method 2 and
% outputs a plot of the coefficient of restitutions sensitivity to varying
% amounts of uncertainty in its inputs assuming a standard bounce time
% vector.
%
% Created by: Johnathan Tucker
% 
% Inputs: figure_number: Figure number for organized plots
%
% Outputs: A plot of the sensitivity of the coefficient of restitution
%          sens_vec_t_prev: A vector of the sensitivity of e to the time of
%                           the previous bounce in seconds
%          sens_vec_tn: A vector of the sensitivity of e to the current
%                       bounce time in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sens_vec_t_prev,sens_vec_tn] = sens_method_2(figure_number)
% Time between bounces in seconds
time_vec = [20,18,16,14,13,12,11,10,9,8,8,7,7,7,6,5]/30;
% Uncertainty in seconds
time_unc_vec = linspace(0,1,15);
% Loop through the length of the time vector
for i = 2:length(time_vec)
    % Calculate the error varying the uncertainty in the nth bounces
    sens_vec_tn(i-1) = sqrt( (1/time_vec(i-1) * time_unc_vec(i-1))^2 + ...
            (-time_vec(i)/(time_vec(i-1)^2) * 0.05)^2);
    % Calculate the error varying the uncertainty in the previous bounce
    sens_vec_t_prev(i-1) = sqrt( (1/time_vec(i-1) * 0.05)^2 + ...
            (-time_vec(i)/(time_vec(i-1)^2) * time_unc_vec(i-1))^2);
end
% Create the figure
figure(figure_number)
% Plot the nth bounce sensitivity
plot(time_unc_vec,sens_vec_tn,'LineWidth',2)
hold on
% Plot the previous bounce sensitivity
plot(time_unc_vec,sens_vec_t_prev,'LineWidth',2)
title('$Sensitivity\:of\:Coefficient\:of\:Restitution\:to\:Method\:2$',...
    'Interpreter','latex',...
    'FontSize',26)
xlabel('$Bounce\:Time\:Uncertainty\:(s)$','Interpreter','latex','FontSize',26)
ylabel('$Coefficient\:of\:Restitution\:Sensitivity$','Interpreter','latex','FontSize',26)
legend('Time "n" ','Time "n-1"')

end