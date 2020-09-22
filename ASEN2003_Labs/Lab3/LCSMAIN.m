%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main driver script for the locomotive crankshaft lab. It
% instantiates variables, calls functions, and creates plots.
%
% Created by: Johnathan Tucker
%
% Inputs: 
%           Sec011_Group15_7.txt = This is the 7V test data obtained from the
%           experimental setup
%           
%           Sec011_Group15_9.txt = This is the 9V test data obtained from the
%           experimental setup
%
%           Sec011_Group15_11.txt = This is the 11V test data obtained from the
%           experimental setup
%
% Outputs:
%           - The mean and standard deviation of the difference between the
%           experimentally obtained values and the values obtained from the
%           model are printed to the terminal for 7V, 9V, and 11V tests.
%
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear;
close all;
clc;
%% Create variables
radius = .077;
distance = .153;
length_exp = .255;
theta_check = linspace(0,2*pi);
V_check = LCSMODEL(radius,distance,length_exp, theta_check,1300);
theta_check = theta_check*(180/pi);
%% Plot the model checking values
figure(1)
plot(theta_check, V_check, 'LineWidth', 2)
title('$Velocity\:Versus\:Angular\:Position\:of\:Model\:Checking\:Data$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(deg)$','Interpreter','latex','FontSize',26)
ylabel('$Velocity(cm/s)$','Interpreter','latex','FontSize',26)
legend('Model Check Data')
%% Call the data processing function
% Get the experimental data from the files
[theta_exp_7V, w_exp_7V, v_exp_7V] = LCSDATA('Sec011_Group15_7.txt');
[theta_exp_9V, w_exp_9V, v_exp_9V] = LCSDATA('Sec011_Group15_9.txt');
[theta_exp_11V, w_exp_11V, v_exp_11V] = LCSDATA('Sec011_Group15_11.txt');

%% Plot the v_exp versus theta_exp for the 7 volt data
figure(2)
plot(theta_exp_7V, v_exp_7V, 'LineWidth', 2)
title('$Velocity\:Versus\:Angular\:Position\:of\:7\:Volt\:Experiment\:Data$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(deg)$','Interpreter','latex','FontSize',26)
ylabel('$Velocity(cm/s)$','Interpreter','latex','FontSize',26)
ylim([-200,250]);
xlim([200,2200]);
legend('Experiment Data')

%% Plot the v_exp versus theta_exp for the 9 volt data
figure(3)
plot(theta_exp_9V, v_exp_9V, 'LineWidth', 2)
title('$Velocity\:Versus\:Angular\:Position\:of\:9\:Volt\:Experiment\:Data$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(deg)$','Interpreter','latex','FontSize',26)
ylabel('$Velocity(cm/s)$','Interpreter','latex','FontSize',26)
ylim([-200,250]);
xlim([200,2200]);
legend('Experiment Data')

%% Plot the v_exp versus theta_exp for the 11 volt data
figure(4)
plot(theta_exp_11V, v_exp_11V, 'LineWidth', 2)
title('$Velocity\:Versus\:Angular\:Position\:of\:11\:Volt\:Experiment\:Data$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(deg)$','Interpreter','latex','FontSize',26)
ylabel('$Velocity(cm/s)$','Interpreter','latex','FontSize',26)
ylim([-200,250]);
xlim([200,2200]);
legend('Experiment Data')

%% Run the model lcs function and plot the results for the data
% 7V first
% Convert to radians
theta_exp_7V = (theta_exp_7V*(pi/180));
w_exp_7V = (w_exp_7V*(pi/180));
% Call the function
[v_model_7V] = LCSMODEL(radius,distance,length_exp,theta_exp_7V,w_exp_7V);
% Convert back to degrees
theta_exp_7V = (theta_exp_7V*(180/pi));
w_exp_7V = (w_exp_7V*(180/pi));

% Now 9V
% Convert to radians
theta_exp_9V = (theta_exp_9V*(pi/180));
w_exp_9V = (w_exp_9V*(pi/180));
% Call the function
[v_model_9V] = LCSMODEL(radius,distance,length_exp,theta_exp_9V,w_exp_9V);
% Convert back to degrees
theta_exp_9V = (theta_exp_9V*(180/pi));
w_exp_9V = (w_exp_9V*(180/pi));

% Now 11V
% Convert to radians
theta_exp_11V = (theta_exp_11V*(pi/180));
w_exp_11V = (w_exp_11V*(pi/180));
% Call the function
[v_model_11V] = LCSMODEL(radius,distance,length_exp,theta_exp_11V,w_exp_11V);
% Convert back to degrees
theta_exp_11V = (theta_exp_11V*(180/pi));
w_exp_11V = (w_exp_11V*(180/pi));

%% Create plots of model velocity and experimental velocity on same graph for 7V
figure(5)
plot(theta_exp_7V, v_exp_7V, 'LineWidth', 2)
hold on
plot(theta_exp_7V, v_model_7V, '--', 'LineWidth', 2)
title('$Comparison\:of\:Velocity\:Versus\:Angular\:Position\:of\:7\:Volt\:Data$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(deg)$','Interpreter','latex','FontSize',26)
ylabel('$Velocity(cm/s)$','Interpreter','latex','FontSize',26)
ylim([-200,300]);
xlim([200,2200]);
legend('Experiment Data', 'Model Data')

%% Create plots of model velocity and experimental velocity on same graph for 9V
figure(6)
plot(theta_exp_9V, v_exp_9V, 'LineWidth', 2)
hold on
plot(theta_exp_9V, v_model_9V , '--', 'LineWidth', 2)
title('$Comparison\:of\:Velocity\:Versus\:Angular\:Position\:of\:9\:Volt\:Data$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(deg)$','Interpreter','latex','FontSize',26)
ylabel('$Velocity(cm/s)$','Interpreter','latex','FontSize',26)
ylim([-200,300]);
xlim([200,2200]);
legend('Experiment Data', 'Model Data')

%% Create plots of model velocity and experimental velocity on same graph for 11V
figure(7)
plot(theta_exp_11V, v_exp_11V, 'LineWidth', 2)
hold on
plot(theta_exp_11V, v_model_11V, '--', 'LineWidth', 2)
title('$Comparison\:of\:Velocity\:Versus\:Angular\:Position\:of\:11\:Volt\:Data$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(deg)$','Interpreter','latex','FontSize',26)
ylim([-200,300]);
xlim([200,2200]);
ylabel('$Velocity(cm/s)$','Interpreter','latex','FontSize',26)
legend('Experiment Data', 'Model Data')
%% Difference between model and experiment
V7_diff =  abs(v_exp_7V - v_model_7V);
V9_diff = abs(v_exp_9V - v_model_9V);
V11_diff = abs(v_exp_11V - v_model_11V);

% Get the mean for each voltage
mean_7V = ones(length(theta_exp_7V),1)*mean(V7_diff);
mean_9V = ones(length(theta_exp_9V),1)*mean(V9_diff);
mean_11V = ones(length(theta_exp_11V),1)*mean(V11_diff);
%% Plot the difference for 7V
figure(8)
plot(theta_exp_7V,abs(V7_diff))
hold on
plot(theta_exp_7V,mean_7V , '--', 'LineWidth', 2)
title('$Difference\:Between\:Modeled\:and\:Experimental\:Velocity:of\:7\:Volt\:Experiment\:Data$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(deg)$','Interpreter','latex','FontSize',26)
ylabel('$Velocity\:Difference(cm/s)$','Interpreter','latex','FontSize',26)
legend('Difference', 'Mean of Difference')
%% Plot the difference for 9V
figure(9)
plot(theta_exp_9V,abs(V9_diff))
hold on
plot(theta_exp_9V, mean_9V, '--', 'LineWidth', 2)
title('$Difference\:Between\:Modeled\:and\:Experimental\:Velocity:of\:9\:Volt\:Experiment\:Data$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(deg)$','Interpreter','latex','FontSize',26)
ylabel('$Velocity\:Difference(cm/s)$','Interpreter','latex','FontSize',26)
legend('Difference', 'Mean of Difference')
%% Plot the difference for 11V
figure(10)
plot(theta_exp_11V,abs(V11_diff))
hold on
plot(theta_exp_11V,mean_11V, '--', 'LineWidth', 2)
title('$Difference\:Between\:Modeled\:and\:Experimental\:Velocity:of\:11\:Volt\:Experiment\:Data$',...
    'Interpreter','latex','FontSize',26)
xlabel('$Angular\:Position\:Theta(deg)$','Interpreter','latex','FontSize',26)
ylabel('$Velocity\:Difference(cm/s)$','Interpreter','latex','FontSize',26)
legend('Difference', 'Mean of Difference')
%% Display the mean and std dev of each difference
fprintf("The mean and std dev of 7V experiment is : %f and %f\n",mean_7V(1,1),std(V7_diff));
fprintf("The mean and std dev of 9V experiment is : %f and %f\n",mean_9V(1,1),std(V9_diff));
fprintf("The mean and std dev of 11V experiment is : %f and %f\n",mean_11V(1,1),std(V11_diff));