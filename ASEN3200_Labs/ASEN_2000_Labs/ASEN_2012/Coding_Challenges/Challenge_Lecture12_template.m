%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE for Lecture 12 - Linear Least-Squares Fit
%
% The purpose of this program is to calculate the equation of the best fit
% line for a data set using linear least-squares fitting.
%
% To complete the challenge, finish the code below to:
% 1) load data from .mat file
% 2) find linear best fit coefficients
% 3) plot the original data along with the best fit line for time t=0 to
%    time t=150s
% 4) add errorbars for fit uncertainty to this plot
%
% Please upload the finished script to Canvas to finish this coding
% challenge.
% 
% STUDENT TEAMMATES
% 1.Johnathan Tucker
% 2.Alex Hill
% 3.
% 4.
%
% CHALLENGE AUTHOR
% Justin Fay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping: Clear the variables and close all open plots
clearvars   % Clears all the variables
close all   % Closes all the plots
clc         % Clears the command line

%% Load data

% Load 'data' structure from .mat file
load('./data.mat')

% Seperate time and y into individual vectors
t = data.t; % [s]
y = data.y; % [deg C]

% Find number of data points
N = numel(t);

%% Find linear best fit coefficients m and b
% Only use polyfit to check other methods
A = ones(length(t),2);
for i = 1:length(A)
    A(i,1) = t(i);
end
x = A\y;
m = x(1,1);
b = x(2,1);
%% Find uncertainty associated with this fit line
%Find sigma y
sum = 0;
for i = 1:length(t)
    temp = (y(i,1)-m*t(i,1)-b)^2;
    sum = sum + temp;
end
N = length(t)-2;
sigma_y = sqrt((1/N)*sum);
%Create W matrix
for i = 1:length(t)
    v(1,i) = (1/(sigma_y^2));
end
W = diag(v);
%Solve for Q and sigma m and sigma b
Q = inv((A')*(W)*(A));
sigma_m = sqrt(Q(1,1));
sigma_b = sqrt(Q(2,2));

%% Plot linear fit line between t=0 and t=150s

tnew = 0: 1 : 150;
line = m*tnew + b;
for i = 1:length(tnew)
    error(1,i) = [tnew(1,i),1]*Q*[tnew(1,i);1];
end

scatter(t,y)
hold on
plot(tnew,line)
hold on
errorbar(tnew,line,error)
title('Given Data With Best Fit Line')
xlabel('Time (s)')
ylabel('Temperature (C)')
