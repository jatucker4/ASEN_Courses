%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL CODING CHALLENGE - PARTIAL PIVOTING/GAUSSIAN ELIMINATION
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1
k_1 = 35;
k_2 = 20.5;
k_3 = 10;
k_4 = 25;
g = 9.81;
m_1 = 3;
m_2 = 4;
m_3 = 1.2;
% Create coefficient matrix
A = [(k_1+k_2+k_3), -k_3, 0 ; -k_3, (k_3 + k_4), -k_4; 0, -k_4, k_4];
b = [m_1*g ; m_2*g ; m_3*g];
u_matrix = A\b;
disp('The u_matrix is:')
disp(u_matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2
% Create A matrix
A_2 = [0, -0.7071, -1, 0, 0, 0.6585, 1, 0, 0, 0, 0, 0, 0, 0;...
    -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
    0, -0.7071, 0, -1, 0, -0.7525, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0.6585, 1, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, -1, 0, -0.7525, 0, 0, 0, 0;...
    0, 0, 0, 0, -1, -0.6585, 0 ,0, 1, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0.7525, 0, 1, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0.7071;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -0.7071;...
    0, 0, 0, 0, 0, 0, 0, 0, -1, -0.6585, 0, 0, 1, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7525, 0, 1, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -0.7071;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7071];
b_2 = [0; 0; 1000; 1800; 0; 800; 0; 5000; 0; 600; 0; 0; 0; 2000];
Forces = A_2\b_2;
largest_force_member = find(Forces == max(Forces));
smallest_force_member = find(Forces == min(Forces));
fprintf('Largest force member is: %0.3g\n', largest_force_member)
fprintf('Smallest force member is: %0.3g\n', smallest_force_member)