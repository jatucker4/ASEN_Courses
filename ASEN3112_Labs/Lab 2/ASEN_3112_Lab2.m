%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3111 - CA2
% 
% Created By: Johnathan Tucker
%
% Collaborators: N/A
%
% The purpose of the script is to contain all of the constants and basic
% computation that will be passed to the plot airfoil flow function as well
% as the plot settings function
%
% Created Date: 2/23/2020
%
% Change Log: 
%           - 2/25/2020: Code up plot airfoil function 
%           - 2/26/2020: Code up plot settings function
%           - 2/27/2020: Finalize formatting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
tic
%Global formatting commands to imporve graphing looks:
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex'); 
%% Begin code block for the first question
% First load the ten pound increment data into a matrix
ten_lb = readmatrix("ASEN_3112_10_lb_increments.csv");

% Process the data for the five pound increment struct.
[F0,F1,F2,F3D,LVDT] = data_process(ten_lb); %[N]

% Create const for the applied loads
applied_load = [0,10,20,30,40,50]*0.453592; %[kg]

% Get the linear fit for each of the datasets
pF0 = polyfit(applied_load,F0,1);
pF1 = polyfit(applied_load,F1,1);
pF2 = polyfit(applied_load,F2,1);
pF3D = polyfit(applied_load,F3D,1);
pLVDT = polyfit(applied_load,LVDT,1);

% Get the coefficient of determination for quantitative linearity
PcoefF0 = corrcoef(applied_load,F0);
PcoefF1 = corrcoef(applied_load,F1);
PcoefF2 = corrcoef(applied_load,F2);
PcoefF3D = corrcoef(applied_load,F3D);
PcoefLVDT = corrcoef(applied_load,LVDT);

% Plot the results and their linear fits
figure
subplot(2,2,1) 
plot(applied_load, F0, 'o','LineWidth',2);
hold on
plot(linspace(applied_load(1),applied_load(end)),...
    (linspace(applied_load(1),applied_load(end)).*pF0(1) + pF0(2)),...
    'LineWidth',2);
title('F0 Sensor Load vs Applied Load');
xlabel('Applied Load [kg]');
ylabel('F0 Sensor[N]');
legend('Experimental Data','Linear Regression','location','nw')

% F1
subplot(2,2,2)
plot(applied_load, F1, 'o','LineWidth',2);
hold on
plot(linspace(applied_load(1),applied_load(end)),...
    (linspace(applied_load(1),applied_load(end)).*pF1(1) + pF1(2)),...
    'LineWidth',2);
title('F1 Sensor Load vs Applied Load');
xlabel('Applied Load [kg]');
ylabel('F1 Sensor[N]');
legend('Experimental Data','Linear Regression','location','nw')

% F2
subplot(2,2,3)
plot(applied_load, F2, 'o','LineWidth',2);
hold on
plot(linspace(applied_load(1),applied_load(end)),...
    (linspace(applied_load(1),applied_load(end)).*pF2(1) + pF2(2)),...
    'LineWidth',2);
title('F2 Sensor Load vs Applied Load');
xlabel('Applied Load [kg]');
ylabel('F2 Sensor[N]');
legend('Experimental Data','Linear Regression','location','nw')

% F3D
subplot(2,2,4)
hold on 
plot(applied_load, F3D, 'o','LineWidth',2);
plot(linspace(applied_load(1),applied_load(end)),...
    (linspace(applied_load(1),applied_load(end)).*pF3D(1) + pF3D(2)),...
    'LineWidth',2);
hold off
title('F3D Sensor Load vs Applied Load');
xlabel('Applied Load [kg]');
ylabel('F3D Sensor[N]');
legend('Experimental Data','Linear Regression','location','nw')

% LVDT
figure
plot(applied_load, LVDT, 'o','LineWidth',2);
hold on 
plot(linspace(applied_load(1),applied_load(end)),...
    (linspace(applied_load(1),applied_load(end)).*pLVDT(1) + pLVDT(2)),...
    'LineWidth',2);
hold off
title('LVDT vs Applied Load');
xlabel('Applied Load [kg]');
ylabel('LVDT Sensor[m]')
legend('Experimental Data','Linear Regression','location','nw')

