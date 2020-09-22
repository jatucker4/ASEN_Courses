function [e_vec, sigma_e,statistics] = e_method_1(height_array,figure_number,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses the first method for calculating the coefficient of
% restitution (Height based)
%
% Created by: Johnathan Tucker
%
% Inputs: unc_height: Uncertainty in the height measurements
%         height_vec: Vector of heights going from initial to nth           
%         height_array: Cell array of heights and their uncertainties
%         figure_number: The figure number of the plot
%
% Outputs: e_vec: Vector of coefficient of restitution valus
%          sigma_e: Uncertainty in coefficient of restitution values
%          statistics: Vector of the mean, median, std dev, and SEOM
%          Also displays a plot of the coefficient of restitution with
%          error bars
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create variables
trial_vec = 1:10;
e_vec_1 = [];
sigma_e_vec = [];
% Loop through the length of the cell array
for j = 1:length(height_array)
    % Account for radius based on if its the ping pong or yellow ball
    if flag == 1
        % Create height vector
        height_vec = height_array{j,1} - 0.94;
    else
        % Create height vector
        height_vec = height_array{j,1} -0.78;
    end
    % Create height uncertainty vector
    unc_height = height_array{j,2};
    % Loop through every height in the vector
    for i = 2:length(height_vec)+1
        % Calculate the coefficient of restitution
        e_vec(i-1) = (height_vec(i-1)/height_vec(1))^(1/(2*(i-1)));
        % Calculate the error in the coefficient of restitution
        sigma_e(i-1) = sqrt(( 1/(2*(i-1)) * (1/height_vec(1)) * ...
            (height_vec(i-1)/height_vec(1))^((1/(2*(i-1)))-1) *...
            unc_height(i-1))^2  + ( 1/(2*(i-1)) * ...
            (-height_vec(i-1)/(height_vec(1)^2)) * ...
            ((height_vec(i-1)/height_vec(1))^((1/(2*(i-1)))-1)) *...
            unc_height(1))^2); 
    end
    % Create a 10 by 6 vector for the 10 trials and 6 measurements in each
    e_vec_1 = [e_vec_1;e_vec];
    % Do the same for the uncertainty
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
%% Create Plots
% If it's the tracker data
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
    hold on
    plot(0:4,statistics(1,4))
    title('$Coefficient\:of\:Restitution\:Method\:One$','Interpreter','latex',...
        'FontSize',26)
    xlabel('$Trial\:Number$','Interpreter','latex','FontSize',26)
    xlim([0 4]);
    xticklabels({'','','1','','2','','3',''});
    ylabel('$Coefficient\:of\:Restitution$','Interpreter','latex','FontSize',26)
    legend('Mean Coefficient of Restitution',...
        'Coefficient of Restitution With Error', 'Standard Error of the Mean')
    hold off
else
    % Create figure number
    figure(figure_number)
    % Create mean vector
    mean_vec = ones(1,12)*mean(e_vec);
    % Plot mean vector
    plot(0:11,mean_vec,'--')
    hold on
    % Plot error bars and scatter plot
    errorbar(trial_vec,e_vec,sigma_e,'o')
    title('$Coefficient\:of\:Restitution\:Method\:One$','Interpreter','latex',...
        'FontSize',26)
    xlabel('$Trial\:Number$','Interpreter','latex','FontSize',26)
    xlim([0 11]);
    xticklabels({'','1','2','3','4','5','6','7','8','9','10',''});
    ylabel('$Coefficient\:of\:Restitution$','Interpreter','latex','FontSize',26)
    legend('Mean Coefficient of Restitution',...
        'Coefficient of Restitution With Error')
    hold off
end
end

