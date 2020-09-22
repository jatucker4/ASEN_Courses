%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 2 - Calculate and plot trajectory of elliptical orbits
%
% The purpose of this program is to calculate the trajectory of various
% objects orbiting Earth and plot them together.
%
% To complete this challenge you will need to complete both this script and
% the sub function ellipticalOrbit.m by doing the following:
% 1) add the elliptical orbit equation to ellipticalOrbit.m
% 2) use ellipticalOrbit.m to calculate the trajectory for multiple orbits
% 3) plot all trajectories in a single figure
% 4) properly label all axes and include a title
% 5) add a legend to label each orbit
%
% If you complete all of these tasks before the class is over, try
% importing the orbital parameters from orbits.csv instead of manually
% entering them.
%
% Please ZIP and upload the ellipticalOrbit function and your team's
% script(s) to Canvas to complete the challenge.
% 
% STUDENT TEAMMATES
% 1. Johnathan Tucker
% 2. Donald Wolfe
% 3. John Wissler
% 4. Dominic Dougherty
%
% CHALLENGE AUTHORS
% Justin Fay 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping: Clear the variables and close all open plots
      % Clears all the variables
      % Closes all the plots
      % Clears the command line
      clear; close all; clc; 
%
%% Initialization: Define the Functions, Load Data, and Set Known Values
% 
% The orbital parameters are specified for you below as seperate vectors
e = [0; 0.5; 0.5; 0.75];
a = [1000; 1000; 1500; 1000];
b = [];
theta = linspace(0,360, 200);
trajectories = [];
% EXTRA: load orbital data from csv file

%% Calculate trajectories for each orbit
%
% Use ellipticalOrbit to calculate a trajectory for each set of paramters.
% These trajectories should include a single full orbit around Earth. If
% you choose to use a discrete set of angles, you should use at least 200.
for i =1:4
    for j = 1:200
        trajectories(j) = ellipticalOrbit(e(i), a(i), theta(j));
    end
    b(i) = a(i)*sqrt(1-e(i)^2);
end


%% Plot results
%
% You should plot all trajectories on the same plot, so they are easy to
% compare. The plots should be in cartesian coordinates with properly
% labeled axes (for this challenge the axes can just be X and Y). A legend
% should be included that clearly shows what orbital parameters (e and a)
% were used for each trajectory. You should also plot a single mark at the
% origin to clearly show the position of Earth.
x = [];
y = [];
polarplot(theta, trajectories)
% for i = 1:4
%     polarplot(theta(i), trajectories(i))
%     hold on
% end

plot(x,y)

