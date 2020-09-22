function [time,theta,omega] = dataProcess(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function process the data from the labview visual interface by
% limiting the angular position values to be between 0.5 and 15 radians
%
% Inputs:
%           filename: Labview file
%
% Outputs: 
%           time: The time the wheel took to travel the angular position
%                 (s)
%           theta: Angluar position of the wheel (radians)
%           omega: Angular velocity of the wheel (radians/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in data
data = load(filename);
theta = data(:,2);
omega = data(:,3);
time = data(:,1);

% cut off theta values
limit1 = 0.5;
limit2 = 15;
spot1 = find(theta > limit1,1);
spot2 = find(theta > limit2,1);
theta = theta(spot1:spot2-1);

omega = omega(spot1:spot2-1);
time = time(spot1:spot2-1);
end