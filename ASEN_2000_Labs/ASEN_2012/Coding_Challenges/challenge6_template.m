%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 6 - Numerical Integration Basics
% 
% This script includes defintions for the trapezoidal method and Simpson's
% 1/3 method. Please read the instructions carefully, and complete the
% three steps outlined below. Bonus questions are included at the end.
%
% Please ZIP and upload your team's script(s) and figures to Canvas to 
% complete the challenge.
% 
% STUDENT TEAMMATES
% 1. Johnathan Tucker
% 2. Alex Hill
% 3.
% 4.
%
% CHALLENGE AUTHORS
% Jelliffe Jackson, Allison Anderson, Torin Clark, John Jackson 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clear all
close all
clc

%% Flow Velocity Measurements in a Pipe

% You have measurements of the laminar flow velocity at different pipe radii:
r = [0.0 0.25 0.5 0.75 1 1.25 1.5 1.75 2.0 2.25]';
v = [38.0 37.6 36.2 33.6 29.7 24.5 17.8 9.6 4.3 0]';

% You want to determine the volumetric flow rate through the pipe.
% The volumetric flow rate through a pipe is analytically found by the
% equation Q = integral(2pi * v(r) * r dr) from 0 to r. However, you've
% taken data on the flow velocity, so you want to perform a numerical
% integration.

figure(1)
plot(v, r, 'k*', 'LineWidth',2);
title('Velocity of Flow in Pipe')
ylabel('r (in)')
xlabel('v (in/s)')
hold on; grid on;

%% (1) Calculate the values of f(r) and make a plot showing f(r) vs r
f_r = 2*pi*v.*r;% Hint: What is the function you're integrating?
figure(2)
plot(f_r,r)
title('f(r) versus r')
xlabel('r (in)')
ylabel('f_r (in^2/s)')
%% (2) Finish defining the function trapezoidal and calculate Q (below)
fprintf('The Trapezoid rule function yields: %f\n',trapezoidal(r,f_r));
% Use the MATLAB function trapz to check your code.
fprintf('The reality check for the trapezoidal rule is: %f\n',trapz(r,f_r));
%% (3) Finish defining the function simpson_onethird and calculate Q (below)
Q_2 = simpson_onethird(r,f_r);
fprintf('The Simpsons Rule functions yields: %f\n', Q_2)
% How can you reality check the result for Simpson's 1/3 method?
simpsons_check = (Q_2 -trapezoidal(r,f_r));
fprintf('The reality check for simpsons rule: %f\n',simpsons_check);

%% NOTE: Functions in a script need to be at the bottom. Don't put any code
% below these functions.
function res = trapezoidal(X, FX)
    sum = 0;
    for i = 2:length(X)
        sum = sum + (X(i) -  X(i-1))*(FX(i) + FX(i-1));
    end
    res = 0.5*sum;    % Hint: you know x and f(x)
end

function res = simpson_onethird(X, FX)
    h = (X(end)-X(1))/length(X);
    evens = 2*sum(FX(2:2:end));
    odds = 4*sum(FX(1:2:end));
    res = (1/3)*h*(FX(1) + evens + odds + FX(end));    % Hint: you know x and f(x)
end

%% Bonus Questions
% In this case, which integration method is better for estimating the
% volumetric flow through the pipe?

% How might you propogate uncertainties in the velocity measurements in the
% numerical integration?
