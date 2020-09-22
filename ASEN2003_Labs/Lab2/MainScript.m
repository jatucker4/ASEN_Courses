%% MainScript.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ASEN 2003 Lab 2
%Group 7
%
%Purpose: This script calls all of the sub functions rquired to execute
%Group 7's code for lab 2. 
%
% Created by: Michael Martinson
%
% Inputs: None        
%
% Outputs:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Recorded Data
clear; clear all; clc;

%Process Data for methods 1, 2, 3:
[iPhoneHeightdata, iPhoneTimeData,StopwatchStoptdata,YellowHeightData] =  ... 
    ProcessData('FullData.xlsx','TimeToStopReadings.xlsx',...
    'Yellow Ball Analysis.xlsx');
 
% Get the coefficient of restitution and plot for method 1
[e_method_1_iphone,e_unc_method_1,statistics_method_1] = e_method_1(iPhoneHeightdata,1,0);
% Get the coefficient of restitution and plot for method 2
[e_method_2_iphone,e_unc_method_2,statistics_method_2] = e_method_2(iPhoneTimeData,2,0);
% Get the coefficient of restitution and plot for method 3
[e_method_3_iphone,e_unc_method_3,statistics_method_3] = e_method_3(StopwatchStoptdata(:,2),...
    iPhoneHeightdata,StopwatchStoptdata(:,1),3,0);
%% Coefficients of Restitution with improved method 1

% Call method 1 on the yellow ball data
[e_method_1_yellow,e_unc_method_1_yellow,statistics_method_1_yellow] = e_method_1(YellowHeightData,4,1);

% Process the improved method data
[ImprovedData, iPhoneTimeData,StopwatchStoptdata,YellowHeightData] =  ... 
    ProcessData('Tests1-3.xlsx','TimeToStopReadings.xlsx',...
    'Yellow Ball Analysis.xlsx');
% Call method 1 on the improved method data
[e_method_1_improved,e_unc_method_1_improved,statistics_method_1_improved] = e_method_1(ImprovedData,5,0);

% Process the Tracker Data
[TrackerData, TrackerTimeData,TrackerStopData,YellowHeightData] =  ... 
    ProcessData('Tracker.xlsx','TimeToStopTracker.xlsx',...
    'Yellow Ball Analysis.xlsx');
% Truncate the tracker data to the necessary trials
TrackerData = TrackerData(1:3,:);
TrackerTimeData = TrackerTimeData(1:3,:);
TrackerStopData = TrackerStopData(1:3,:);
% Remove the extra zero from every trials vector
for i = 1:length(TrackerData)
    temp_vec = TrackerData{i,1};
    temp_vec_2 = TrackerData{i,2};
    temp_vec(6,:) = [];
    temp_vec_2(6,:) = [];
    TrackerData{i,1} = temp_vec;
    TrackerData{i,2} = temp_vec_2;
end

% Call method 1 on the tracker data
[e_method_1_tracker,e_unc_method_1_tracker,statistics_method_1_tracker] = e_method_1(TrackerData,6,2);
% Call method 2 on the tracker data
[e_method_2_tracker,e_unc_method_2_tracker,statistics_method_2_tracker] = e_method_2(TrackerTimeData,7,2);
% Call method 3 on the tracker data
[e_method_3_tracker,e_unc_method_3_tracker,statistics_method_3_tracker] =...
    e_method_3(TrackerStopData(:,2),TrackerData,TrackerStopData(:,1),8,2);
%% Make sensitivity function calls
% Get the sensitivity of method 1
[sens_e_method_1_h0,sens_e_method_1_hn] = sens_method_1(9);
% Get the sensitivity of method 2
[sens_e_method_2_t_prev,sens_e_method_2_tn] = sens_method_2(10);
% Get the sensitivity of method 3
[sens_e_method_3_t_total,sens_e_method_3_h0] = sens_method_3(11);