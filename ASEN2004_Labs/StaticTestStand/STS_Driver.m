%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver script for the static test stand lab. It handles all
% function calls and plotting
%
% Written by: Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear;
close all;
clc;
%% Navigate to directory and get the data
% This is a structure with information on every file in the directory
data_directory = dir('STSData');
for i = 1:length(data_directory)
    % Only pull the data from the ones with Group in their name
    if contains(data_directory(i,1).name,'Group')
        % Get the path to the data file
        path = string(data_directory(i,1).folder);
        % Get the file name
        filename = string(data_directory(i,1).name);
        % Obligatory slash
        filename = strcat(path,"\",filename);
        % Put it all together
        data{i,:} = sts_dataprocess(filename);
    end
end
% Make sure everything is closed to keep matlab happy
fclose('all');
% Remove hidden folders that matlab found in the directory
data = data(3:end);
%% Pull the thrust data from each test and plot it
for i = 1:length(data)
    thrust = data{i};
    thrust = thrust{6};
    time_vec = data{i};
    time_vec = time_vec{8};
    total_time(i) = time_vec(end) - time_vec(1);
    max_thrust(i) = max(thrust);
    hold on
    plot(time_vec,thrust,'b');
end
hold on
% hline(mean(max_thrust),'-','b');
xlabel("$Time\:(s)$",'FontSize',20,'Interpreter',...
    'Latex');
ylabel("$Thrust\:(N)$",'FontSize',20,'Interpreter',...
    'Latex');
title("$Thrust\:Data\:From\:Static\:Test\:Stand\:Trials$",...
    'FontSize',20,'Interpreter','Latex');
legend("Thrust Data");

% Create peak thrust histogram
figure(2)
histogram(max_thrust);
xlabel("$Peak\:Thrust\:(N)$",'FontSize',20,'Interpreter',...
    'Latex');
ylabel("$Number\:of\:Trials$",'FontSize',20,'Interpreter',...
    'Latex');
title("$Histogram\:of\:Peak\:Thrust\:Values$",...
    'FontSize',20,'Interpreter','Latex');

fprintf("The average thrust value is: %f Newtons\n",mean(max_thrust));
fprintf("The standard deviation of thrust is: %f Newtons\n",std(max_thrust));

%% Get ISP Data
for i = 1:length(data)
    isp_temp = data{i};
    isp(i) = isp_temp{7};
end
figure(3)
histogram(isp);
xlabel("$ISP\:(s)$",'FontSize',20,'Interpreter',...
    'Latex');
ylabel("$Number\:of\:Trials$",'FontSize',20,'Interpreter',...
    'Latex');
title("$Histogram\:of\:ISP\:Values$",...
    'FontSize',20,'Interpreter','Latex');

fprintf("Average ISP is: %f seconds\n",mean(isp));
fprintf("The standard deviation is: %f seconds\n",std(isp));
fprintf("The standard error of the mean for ISP is: %f seconds\n",...
    std(isp)/sqrt(length(isp)));

%% Plot the total time histogram
figure(4)
histogram(total_time);
xlabel("$Total\:Thrusting\:Time\:(s)$",'FontSize',20,'Interpreter',...
    'Latex');
ylabel("$Number\:of\:Trials$",'FontSize',20,'Interpreter',...
    'Latex');
title("$Histogram\:of\:Total\:Thrust\:Time\:Values$",...
    'FontSize',20,'Interpreter','Latex');

fprintf("The average total thrust time value is: %f seconds\n",...
    mean(total_time));
fprintf("The standard deviation of total thrust time is: %f seconds\n",...
    std(total_time));

%% Plot of SEM vs N and Confidence Interval Stuff
N = 1:1:1000;
SEM_vec = std(isp)./sqrt(N);
SEM_vec_95 = SEM_vec.*1.96;
SEM_vec_975 = SEM_vec.*2.24;
SEM_vec_99 = SEM_vec.*2.58;

% Get the number of trials necessary for each respective CI in a 0.1s mean
% ISP
N_95_1 = find(SEM_vec_95 < 0.1, 1, 'first');
N_975_1 = find(SEM_vec_975 < 0.1, 1, 'first');
N_99_1 = find(SEM_vec_99 < 0.1, 1, 'first');

% Get the number of trials necessary for each respective CI in a 0.01s mean
% ISP
N_95_01 = find(SEM_vec_95 < 0.01, 1, 'first');
N_975_01 = find(SEM_vec_975 < 0.01, 1, 'first');
N_99_01 = find(SEM_vec_99 < 0.01, 1, 'first');

% Plot the confidence interval vectors
figure(5)
plot(N,SEM_vec,'LineWidth',2);
hold on
plot(N,SEM_vec_95,'LineWidth',2);
hold on
plot(N,SEM_vec_975,'LineWidth',2);
hold on
plot(N,SEM_vec_99,'LineWidth',2);
xlabel("$SEM\:(s)$",'FontSize',20,'Interpreter',...
    'Latex');
ylabel("$Number\:of\:Trials$",'FontSize',20,'Interpreter',...
    'Latex');
title("$SEM\:vs\:Number\:of\:Trials$",'FontSize',20,'Interpreter',...
    'Latex');
legend("SEM vs N(No CI)","SEM vs N(95% CI)","SEM vs N(97.5% CI)",...
    "SEM vs N(99% CI)")
% Lots of print statements
fprintf("For a CI of 95%% in a mean of 0.1s it would require %d trials\n",...
N_95_1);
fprintf("For a CI of 97.5%% in a mean of 0.1s it would require %d trials\n",...
N_975_1);
fprintf("For a CI of 99%% in a mean of 0.1s it would require %d trials\n",...
N_99_1);
fprintf("For a CI of 95%% in a mean of 0.01s it would require %d trials\n",...
N_95_01);
fprintf("For a CI of 97.5%% in a mean of 0.01s it would require %d trials\n",...
N_975_01);
fprintf("For a CI of 99%% in a mean of 0.01s it would require %d trials\n",...
N_99_01);