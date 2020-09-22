function [theta,omega] = dataProcess(filename)
% Read in data
data = load(filename);
theta = data(:,2);
omega = data(:,3);

% cut off theta values
limit1 = 0.5;
limit2 = 15;
spot1 = find(theta > limit1,1);
spot2 = find(theta > limit2,1);
theta = theta(spot1:spot2-1);

omega = omega(spot1:spot2-1);
end