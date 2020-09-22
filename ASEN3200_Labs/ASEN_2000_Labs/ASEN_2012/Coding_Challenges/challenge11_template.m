%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 11 - Zero Finding with MATLAB
% 
% This challenge is an exercise in finding the zeros of functions in
% MATLAB. There are several functions you can use to do this - use the
% `help` function to learn more. You will find the zeros of the van der
% Waals equation of a gas state for CO2. 
%
% Please ZIP and upload your team's script(s) and figures to Canvas to 
% complete the challenge.
%
% STUDENT TEAMMATES
% 1.
% 2.
% 3.
% 4.
%
% CHALLENGE AUTHORS
% Allison Anderson, John Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% Given:
% The equation and constants for CO2. a and b are specific to CO2.
a = 3.59; % atm*L^2/mol^2
b = 0.0427; % l/mol
R = 0.082; % L*atm/K*mol
P = 10; % atm
T = 300; % K

vdw = @(V) (P-a./V.^2).*(V-b) - R*T; % The anonymous function you will use.

% Volumes over which to plot/investigate for roots
% 1. Create a volume vector linearly spaced between -1 and 3
volume = linspace(-1,3);
% 2. Plot the van der Waals function vs. volume. Change the scale of the
% plot to get a better look around the roots.
van = vdw(volume);
plot(volume, van)
hold on
% 3. From looking at the plot, how many roots are there, and what are they
% approximately?
% 3 roots
% 4. Find the roots of the polynomial using the `roots` function.
coefficients = [P, -(P*b + R*T), -a, a*b];
zeros = roots(coefficients);
% 5. Find the roots of the polynomial using the `fzero` function.
x0 = [-0.15, 0.05, 2.6];
zeros_3 = fzero(vdw,-0.15);
zero_3 = fzero(vdw,0.05);
zero_4 = fzero(vdw, 2.6);
% 6. Find the roots of the polynomial using the `fsolve` function.

zeros_2 = fsolve(vdw,x0);
% 7. Plot the roots that you calculated on the plot created in 2. Just
% choose one set of the roots (4, 5 OR 6)
vline(zeros_2)
% 8. Calculate the % error between the roots using the different methods.
% Is it an acceptable error?

% Bonus question: Go to this website and choose another gas to test: 
% http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/waal.html
% HINT: Make sure you check your units!

% Bonus question: change the parameters for the CO2 to see how the
% locations of the zeros change. Are there any parameter values that result
% in no zeros for the function?











