function [e_vec, sigma_e,statistics] = e_method_2(time_array,figure_number,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses the second method for calculating the coefficient of
% restitution (Time between bounces)
%
% Created by: Johnathan Tucker
%
% Inputs: unc_time: Uncertainty in the time measurements
%         time_vec: Vector of time between bounce          
%         time_array: Cell array of time trials and their uncertainties
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
e_vec_1 = [];
sigma_e_vec = [];
% Loop through the array
for j = 1:length(time_array)
    % For every trial create a time and uncertainty vector
    time_vec = time_array{j,1};
    unc_time = time_array{j,2};
    % Loop through the time and uncertainty vectors
    for i = 2:length(time_vec)
        % Calculate the coefficient of restitution
        e_vec(i-1) = time_vec(i)/time_vec(i-1);
        % Calculate the uncertainty
        sigma_e(i-1) = sqrt( (1/time_vec(i-1) * unc_time(i))^2 + ...
            (-time_vec(i)/(time_vec(i-1)^2) * unc_time(i-1))^2);
    end
    % Append the calculated values to a vector
    e_vec_1 = [e_vec_1;e_vec];
    sigma_e_vec = [sigma_e_vec;sigma_e];
end
%% Calculate the weighted average
[m,n] = size(sigma_e_vec);
% For every row
for i = 1:m
    % And every column in the row
    for j = 1:n
       % Caclulate the uncertainty at that row/col position
       w_i(j) = 1/(sigma_e_vec(i,j)^2);    
    end
    % Get the weighted average of the coefficient of restitution in the row
    e_vec(i) = sum(w_i.*e_vec_1(i,:))/sum(w_i);
    % Get the uncertainty for that row
    sigma_e(i) = 1/sqrt(sum(w_i));
end
%% Calculate and output statistics for the e_vec
statistics(1,1) = mean(e_vec);
statistics(1,2) = median(e_vec);
statistics(1,3) = std(e_vec);
statistics(1,4) = std(e_vec)/sqrt(length(e_vec));
%% Create plots for outputting
if flag == 2
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
    title('$Coefficient\:of\:Restitution\:Method\:Two\:Tracker\:Data$',...
        'Interpreter','latex','FontSize',26)
    xlabel('$Trial\:Number$','Interpreter','latex','FontSize',26)
    xlim([0 4]);
    xticklabels({'','','1','','2','','3',''});
    ylabel('$Coefficient\:of\:Restitution$','Interpreter','latex','FontSize',26)
    legend('Mean Coefficient of Restitution',...
        'Coefficient of Restitution With Error')
    hold off
else 
    % Create the figure
    figure(figure_number)
    % Create the mean vector
    mean_vec = ones(1,12)*mean(e_vec);
    % Plot the mean vector
    plot(0:11,mean_vec,'--')
    hold on
    % Plot the errorbars and scatter plot
    errorbar(trial_vec,e_vec,sigma_e,'o')
    title('$Coefficient\:of\:Restitution\:Method\:Two$','Interpreter','latex',...
        'FontSize',26)
    xlabel('$Trial\:Number$','Interpreter','latex','FontSize',26)
    xlim([0 11]);
    xticklabels({'','1','2','3','4','5','6','7','8','9','10',''});
    ylabel('$Coefficient\:of\:Restitution$','Interpreter','latex','FontSize',26)
    legend('Mean Coefficient of Restitution',...
        'Coefficient of Restitution With Error')
end
end