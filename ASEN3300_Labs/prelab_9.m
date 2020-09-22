%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3300 - PreLab 9
% 
% Created By: Johnathan Tucker
%
% Collaborators: N/A
%
% The purpose of the script is to create the truth table for question 1 to
% validate my simplest for siren function
%
% Created Date: 3/18/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
%% Create possible state vectors
A = [0; 0; 0; 0; 1; 1; 1; 1];
O = [0; 0; 1; 1; 1; 1; 0; 0];
F = [0; 1; 0; 1; 1; 0; 1; 0];

%% Get siren truth vector
S = siren(A,O,F);
Truth = [A,O,F,S'];
fprintf("     A     O     F     S\n")
disp(Truth)



