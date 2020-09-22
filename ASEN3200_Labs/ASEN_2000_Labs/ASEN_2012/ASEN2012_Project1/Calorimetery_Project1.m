%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALORIMETERY PROJECT 1 - Estimate the specific heat of an unknown sample
%
% The purpose of this project is to use linear regression to estimate
% temperature values and then to propogate the error using a method of our
% choosing, in this case the general method.
% 
% This code has no inputs or outputs, however it does print the specific 
% heat and uncertainty of this measurement to the command window. 
% 
% Date Created: 10/20/2018
% Date Modified: N/A 
% That's right. I write perfect code the first time
%
% Group Members:
%1) Johnathan Tucker 108602275
%2) Donald Wolfe 108683692
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% General Housekeeping
clear;
close all;
clc;
%% Declare variables
calorimeter_spec_heat = 0.214;
sample_mass = 114.400;
calorimeter_mass = 318.3;
%Pull data from text file
file_data = read_sample_data('Sample_B.txt');
%Split it into variables
time = table2array(file_data(:,2));
water_temp = table2array(file_data(:,3));
calorimeter_temp1 = table2array(file_data(:,4));
calorimeter_temp2 = table2array(file_data(:,5));
%Using those get the temperatures of the sample and calorimeter
sample_initial_temp = mean(water_temp);
%Uncertainty of our initial sample measurement
sample_initial_temp_uncertainty = std(water_temp);

%% Lets Create our Least Squares to Find our Temperatures
% Average the thermocouple values together to get an average temperature
% vector
for i = 1:length(time)
calorimeter_mean_temp(i,1) = mean([calorimeter_temp1(i,1),calorimeter_temp2(i,1)]);
end
% Eyeball the time the sample was inserted
index_inserted = find(time == 218.293);
%Use this time to plot a least squares fit and find the initial temperature
p1 = polyfit(time(1:218,1), calorimeter_mean_temp(1:218,1),1);
f1 = polyval(p1,time(1:220,1));
calorimeter_initial_temp2 = f1(end,1);
%Find the index when the maximum temperature is reached
max_index = find(calorimeter_mean_temp == max(calorimeter_mean_temp));
%Create a least squares the extends from this time to the end
p2 = polyfit(time(max_index:end,1),calorimeter_mean_temp(max_index:end,1),1);
%Extrapolate it back to the time the sample was inserted
f2 = p2(1,1)*(time(218:end,1)) + p2(1,2);
%Create a least squares fit for the time the temperature was rising to find
%our average temp
p3 = polyfit(time(230:275,1), calorimeter_mean_temp(230:275,1),1);
f3 = p3(1,1).*time(200:330,1) + p3(1,2);
%Create a variable for our temp high
calorimeter_final_temp2 = f2(1,1);
%Find the average temperature
avg_temp = (calorimeter_final_temp2+calorimeter_initial_temp2) / 2;
%Find the index of the average temperature in the least squares vector
time_avg = find(f3 >= avg_temp);
%Find the time when the average temp happened
time_avg = time_avg(1,1) + 200;
%Find a better approximation for the second temperature using this average
%time
calorimeter_final_temp2 = p2(1,1)*time_avg + p2(1,2);
%set that equal to the sample's final temperature
sample_final_temp2 = calorimeter_final_temp2;
%% Plots
%Plot of the calorimeter temp over time
plot(time,calorimeter_mean_temp)
hold on
%Plot of the first least squares
plot(f1)
hold on
%Plot of the second least squares
plot(time(200:330,1),f3)
hold on
%Plot of the third least squares
plot(time(218:end,1),f2)
hold on
%Plot a vertical line at the time when Temp average occurs
vline(time_avg, 'r:',  'Time of T_a_v_g')
hold on
%Plot a horizontal line at the initial temp
hline(calorimeter_initial_temp2, 'g:', 'T_o')
hold on
%Plot a horizontal line at the second temperature
hline(calorimeter_final_temp2, 'b:', 'T_2')
%Include proper labels,title, and legend
xlabel('Time (s)')
ylabel('Temperature (C)')
title('Temperature Change of a Calorimeter System with Time')
legend('Calorimeter Temperature', 'Least Squares Fit 1', 'Least Squares Fit 2', 'Least Squares Fit 3')
%% Finding the Uncertainty in the first temperature
%We'll use the matrix methods to find the delta y new that is associated
%with the first linear regression, where the first linear regression is
%defined as the one at the bottom of the graph.

%Create our A matrix for the first fit
A_1 = ones(220,2);
A_1(:,1) = time(1:220,1);
%Calculate sigma y for the weigthed matrix
sum = 0;
for i = 1:220
    sum = (calorimeter_mean_temp(i,1) - p1(1,1)*time(i,1) - p1(1,2))^2 + sum;
    t_new1(1,i) = time(i,1);
    t_new1(2,i) = 1;
end
%Calculate the sigma y value
sigma_y_fitone = sqrt((1/(220))*sum);
%Create the weighted matrix for the first fit
W_1 = (1/sigma_y_fitone)*eye(220);
%Create the Q matrix for the first fit
Q_1 = inv((A_1')*W_1*A_1);
%Create the delta y new for the first fit
for j = 1:220
    delta_y_one(j,1) = A_1(j,:)*Q_1*t_new1(:,j);
end
%Therefor the uncertainty in our T_1 measurement is:
uncertainty_t1 = sqrt(delta_y_one(end,1));
%% Finding the Uncertainty in the Temperature 2 measurement
%We'll use the matrix methods to find the delta y new that is associated
%with the second linear regression, where the second linear regression is
%defined as the one at the top of the graph.

%Create our A matrix for the first fit
A_2 = ones(length(time(218:end,1)),2);
A_2(:,1) = time(218:end,1);
%Calculate sigma y for the weigthed matrix
sum = 0;
for i = 1:length(time(218:end,1))
    sum = (calorimeter_mean_temp(i,1) - p2(1,1)*time(i,1) - p2(1,2))^2 + sum;
    t_new2(1,i) = time(i,1);
    t_new2(2,i) = 1;
end
%Calculate sigma y for the second linear regression
sigma_y_fittwo = sqrt((1/(666))*sum);
%Create the weighted matrix
W_2 = (1/sigma_y_fittwo)*eye(666);
%Create the Q matrix
Q_2 = inv((A_2')*W_2*A_2);
%Create the delta y new matrix for this second fit
for j = 1:length(time(218:end,1))
    delta_y_two(j,1) = A_2(j,:)*Q_2*t_new2(:,j);
end
%Therefor the uncertainty in our T_2 measurement is at time average:
uncertainty_t2 = sqrt(delta_y_two(time_avg,1));
%% Now lets calculate the specific heat using these measurements.
%This is in cal/(g*K)
sample_spec_heat2 = (calorimeter_mass*calorimeter_spec_heat*(calorimeter_final_temp2 - calorimeter_initial_temp2))...
    /(sample_mass*(sample_initial_temp - sample_final_temp2));
%Now we have to convert to our sample comparison units of J/(g*K)
sample_spec_heat2 = sample_spec_heat2*(4.2);
%Therefore I think our sample is Copper
%Lets find percent error between the two
percent_error2 = 100 * abs(sample_spec_heat2 - 0.261)/0.261;
%% Calculating the Uncertainty in our Specific Heat measurement
% In order to find the uncertainty in our specific heat we'll use the
% general method

% Create a function for the partial derivative of sample specific heat with
% respect to T naut
partial_cs_To = (calorimeter_mass*calorimeter_spec_heat*(calorimeter_final_temp2 - 1))...
    /(sample_mass*(sample_initial_temp - sample_final_temp2));
% Create a function for the partial derivative of sample specific heat with
% respect to T 1
partial_cs_T1 = (-calorimeter_mass*calorimeter_spec_heat*(calorimeter_final_temp2 - calorimeter_initial_temp2))...
    /(sample_mass*((sample_initial_temp - sample_final_temp2)^2));
% Create a function for the partial derivative of sample specific heat with
% respect to T 2
partial_cs_T2 = (calorimeter_mass*calorimeter_spec_heat/sample_mass)*...
    (((sample_initial_temp - calorimeter_final_temp2)*(1-calorimeter_initial_temp2) - (calorimeter_final_temp2 - calorimeter_initial_temp2)*(sample_initial_temp - 1))...
    /((sample_initial_temp - sample_final_temp2)^2));
% Create a function for the uncertainty in the specific heat
uncertainty_cs = sqrt((partial_cs_To*sample_initial_temp_uncertainty)^2 ...
    + (partial_cs_T1*uncertainty_t1)^2 + (partial_cs_T2*uncertainty_t2)^2);
fprintf('\nThe specific heat for the sample is %f and the uncertainty is %f\n', sample_spec_heat2, uncertainty_cs);