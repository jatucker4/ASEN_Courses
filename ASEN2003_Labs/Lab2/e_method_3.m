function [e_vec, sigma_e,statistics] = e_method_3(total_time_vec, height_array,...
    total_time_unc,figure_number,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses the third method for calculating the coefficient of
% restitution (Height and Total time)
%
% Created by: Johnathan Tucker
%
% Inputs: total_time_unc: Uncertainty in the total time measurements
%         total_time_vec: Vector of total times from drop to settling           
%         height_array: Cell array of heights and their uncertainties
%         figure_number: The figure number of the plot
%
% Outputs: e_vec: Vector of coefficient of restitution valus
%          sigma_e: Uncertainty in coefficient of restitution values
%          statistics: Vector of the mean, median, std dev, and SEOM
%          A plot is created of the coefficient of restitution with error
%          bars versus the trial number
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create necessary variables
trial_vec = 1:10;
g = 386.09;
% Loop through the array
for j = 1:length(height_array)
    % For every trial create height and uncertainty vectors
    height_vec = height_array{j,1};
    unc_height_vec = height_array{j,2};
    % For every trial create initial height and uncertainty values
    initial_height_vec(j) = height_vec(1,1);
    initial_height_unc(j) = unc_height_vec(1,1);
end
% Loop through the created vectors
for i = 1:length(total_time_vec)
    % Calculate the coefficient of restitution
    e_vec(i) = (total_time_vec(i) - sqrt((2 * initial_height_vec(i))/g))/...
        (total_time_vec(i) + sqrt((2 * initial_height_vec(i))/g));
    % Calculate the uncertainty
    sigma_e(i) = sqrt( (((total_time_vec(i) + sqrt((2 * ...
        initial_height_vec(i))/g)) - ((total_time_vec(i) - sqrt((2 *...
        initial_height_vec(i))/g))))/(((total_time_vec(i) + sqrt((2 *...
        initial_height_vec(i))/g)))^2) * total_time_unc(i))^2 + ...
        (((total_time_vec(i) + sqrt((2 *initial_height_vec(i))/g))*(1/g)*...
        ((2*initial_height_vec(i)/g)^(-1/2)) + (total_time_vec(i) - ...
        sqrt((2 *initial_height_vec(i))/g))*(1/g)*...
        ((2*initial_height_vec(i)/g)^(-1/2)))/(((total_time_vec(i) + sqrt((2 *...
        initial_height_vec(i))/g)))^2) * initial_height_unc(i))^2);
end
%% Calculate and output statistics for the e_vec
statistics(1,1) = mean(e_vec);
statistics(1,2) = median(e_vec);
statistics(1,3) = std(e_vec);
statistics(1,4) = std(e_vec)/sqrt(length(e_vec));
%% Create plots for outputting
if flag ==2
    % Delete unnecessary data
    e_vec(:,4:end) = [];
    sigma_e(:,4:end) = [];
    % Create the figure
    figure(figure_number)
    % Create a mean vector
    mean_vec = ones(1,5)*mean(e_vec);
    % Plot the mean vector
    plot(0:4,mean_vec,'--')
    hold on
    % Plot the error bars and scatter plot
    errorbar(1:3,e_vec,sigma_e,'o')
    title('$Coefficient\:of\:Restitution\:Method\:Three\:Tracker\:Data$',...
        'Interpreter','latex','FontSize',26)
    xlabel('$Trial\:Number$','Interpreter','latex','FontSize',26)
    xlim([0 4]);
    xticklabels({'','','1','','2','','3',''});
    ylabel('$Coefficient\:of\:Restitution$','Interpreter','latex','FontSize',26)
    legend('Mean Coefficient of Restitution',...
        'Coefficient of Restitution With Error')
    hold off
else
    % Create the figure number
    figure(figure_number)
    % Create the mean vector
    mean_vec = ones(1,12)*mean(e_vec);
    % Plot the mean vector
    plot(0:11,mean_vec,'--')
    hold on
    % Plot the errorbars with scatter plot
    errorbar(trial_vec,e_vec,sigma_e,'o')
    title('$Coefficient\:of\:Restitution\:Method\:Three$','Interpreter','latex',...
        'FontSize',26)
    xlabel('$Trial\:Number$','Interpreter','latex','FontSize',26)
    xlim([0 11]);
    xticklabels({'','1','2','3','4','5','6','7','8','9','10',''});
    ylabel('$Coefficient\:of\:Restitution$','Interpreter','latex','FontSize',26)
    legend('Mean Coefficient of Restitution',...
        'Coefficient of Restitution With Error')
end
end