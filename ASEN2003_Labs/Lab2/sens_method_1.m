%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When called this function conducts sensitivity analysis on method 1 and
% outputs a plot of the coefficient of restitutions sensitivity to varying
% amounts of uncertainty in its inputs assuming a perfect initial height
% of 20 inches and a standard bounce height vector.
%
% Created by: Johnathan Tucker
% 
% Inputs: figure_number: Figure number for organized plots
%
% Outputs: A plot of the sensitivity of the coefficient of restitution
%          sens_vec_h0: A vector of the sensitivity of e to initial height
%                       in inches
%          sens_vec_hn: A vector of the sensitivity of e to bounce height
%                       in inches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma_e_h0,sigma_e_hn] = sens_method_1(figure_number)
% Create a vector of the different uncertainties in inches
unc_vec = linspace(0,3,14);
% Create the initial height values in inches
initial_height = 28;
% Create the height vector in inches
height_vec = [30,23,14,11,10,8,7,5.5,4.5,3.5,3,2.5,2,1];

for i = 1:length(unc_vec)
    % Calculate the error in e by only varying the uncertainty in the nth
    % bounce height
    sigma_e_hn(i) = sqrt(( 1/(2*(i)) * (1/initial_height) * ...
            (height_vec(i)/initial_height)^((1/(2*(i)))-1) *...
            unc_vec(i))^2  + ( 1/(2*(i)) * ...
            (-height_vec(i)/(initial_height^2)) * ...
            ((height_vec(i)/initial_height)^((1/(2*(i)))-1)) *...
            0.3)^2);
    % Calculate the error in e by only varying the uncertainty in the
    % initial bounce height    
    sigma_e_h0(i) = sqrt(( 1/(2*(i)) * (1/initial_height) * ...
            (height_vec(i)/initial_height)^((1/(2*(i)))-1) *...
            0.3)^2  + ( 1/(2*(i)) * ...
            (-height_vec(i)/(initial_height^2)) * ...
            ((height_vec(i)/initial_height)^((1/(2*(i)))-1)) *...
            unc_vec(i))^2);
end
% Create the figure 
figure(figure_number)
% Plot the nth bounce sensitivity
plot(unc_vec,sigma_e_hn,'LineWidth',2)
hold on
% Plot the initial bounce sensitivity
plot(unc_vec,sigma_e_h0 ,'LineWidth',2)
title('$Sensitivity\:of\:Coefficient\:of\:Restitution\:to\:Method\:1$',...
    'Interpreter','latex',...
    'FontSize',26)
xlabel('$Bounce\:Height\:Uncertainty(in)$','Interpreter','latex','FontSize',26)
ylabel('$Coefficient\:of\:Restitution\:Sensitivity$','Interpreter','latex','FontSize',26)
legend('Bounce Height','Initial Bounce Height')
end