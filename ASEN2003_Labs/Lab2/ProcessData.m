%% ProcessData.m

%ASEN 2003 Lab 2
%Group 7
%
%Purpose:
%
% Created by: Michael Martinson
%
% Inputs:            
%
% Outputs:
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [iPhoneHeightdata, iPhoneTimeData,StopwatchStoptdata,YellHeightdata] = ...
    ProcessData(Filename1,Filename2,Filename3)

%Process methods 1 and 2 data:
[num1,txt1,raw1] = xlsread(Filename1); 


%% Process Height data: 
iPhoneHeightdata = cell(10,2); 

HeightDatarounded = round(num1(:,1), 1); 

%Number of data sets: 
numsets = 10; 

for i = 0:numsets-1
iPhoneHeightdata{i+1,1} = HeightDatarounded(1+7*(i):6+7*(i),1); 
end

%Add height uncertainty:
HuncVec = ones(6,1)*0.3; 

for i = 1:numsets
iPhoneHeightdata{i,2} = HuncVec;
end 


%% Process time data 
iPhoneTimeData = cell(10,2); 


TimeDatarounded = round(num1(:,2), 2); 

for i = 0:numsets-1
iPhoneTimeData{i+1,1} = TimeDatarounded(1+7*(i):5+7*(i),1); 
end

%Add time uncertainty:
TuncVec = ones(5,1)*0.05; 

for i = 1:numsets
iPhoneTimeData{i,2} = TuncVec;
end 

%% 
%Process method 3 data:
[num2,txt2,raw2] = xlsread(Filename2);

%Uncertainty vector for stopwatch method:
Ustopwatchvec = 0.5*ones(10,1); 
stopwatchdatarounded = round(num2,1); 

StopwatchStoptdata = zeros(10,2); 
StopwatchStoptdata(:,2) = stopwatchdatarounded; 
StopwatchStoptdata(:,1) = Ustopwatchvec; 

%% Process Yellow Ball Data:
[num3,txt3,raw3] = xlsread(Filename3);

YellHeightdata = cell(10,2); 

yelldatarounded = round(num3,1); 

%Number of data sets: 
numsets = 10; 

for i = 0:numsets-1
YellHeightdata{i+1,1} = yelldatarounded(1+7*(i):6+7*(i),1); 
end

%Add height uncertainty:
HuncVec = ones(6,1)*0.3;

for i = 1:numsets
YellHeightdata{i,2} = HuncVec;
end 


%Remove erroneus last data point from the last test:
YellHeightdata{10,1}(end) = []; 
YellHeightdata{10,2}(end) = []; 




end