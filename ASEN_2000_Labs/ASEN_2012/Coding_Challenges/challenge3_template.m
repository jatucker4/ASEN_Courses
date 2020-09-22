%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 3 - Calculate mean and standard deviation and plot error
%                    bars from two data sets
%
% The purpose of this program is to assess whether two data sets are the
% same.
%
% To complete the challenge, finish the code below to:
% 1) load csv file
% 2) calculate mean and standard deviation of data set
% 3) plot mean with standard deviation error bars
% 4) set axes to set data points in center of plot
% 5) correctly label axes, title, and legend
%
% Please ZIP and upload the 'func.txt', your team's script(s), and the
% results file to Canvas to complete the challenge.
% 
% STUDENT TEAMMATES
% 1. Johnathan Tucker
% 2. Isaac Goldner
% 3. Travis Griffin
% 4.
%
% CHALLENGE AUTHOR
% Melinda Zavala
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping: Clear the variables and close all open plots
   % Clears all the variables
   % Closes all the plots
   % Clears the command line
   clear; close all; clc;
%% Initialization: Define the Functions, Load Data, and Set Known Values
% load data using csvread (data file format). Start at row 2 and read all
% columns.
data = csvread('subject_data.csv', 1, 0);
subject_1 = data(:,1);
subject_2 = data(:,2);
%% Calculate mean and standard deviations for each subject
%
% calculate mean
subject_1_mean = mean(subject_1);
subject_2_mean = mean(subject_2);
% calculate standard deviation
subject_1_standev = standev(subject_1, subject_1_mean);
subject_2_standev = standev(subject_2, subject_2_mean);

%% Plotting: Plot subject data with error bars
%
% plot mean as a point, and standard deviation error bars
% use 'errobar' function to create bars
names = {'subject_1','subject_2'};
errorbar(1,subject_1_mean, subject_1_standev,'o');
hold on;
errorbar(2,subject_2_mean, subject_2_standev,'o');
set(gca,'xtick',[1:2],'xticklabel', names);
% set x-axes limits so that points show up in center of plot
xlim([0,3])
% label axes, title, and legend
ylabel('Mean of the Subject Data');
xlabel('Subjects');
legend('subject_1', 'subject_2');
title('Subject Mean and Standard Deviation');
% save file as '.png' type. Use 'saveas' function
saveas(gcf, 'subject_mean_stdev_plot.png');
