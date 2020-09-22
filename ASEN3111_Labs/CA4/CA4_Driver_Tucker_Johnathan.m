%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3111 - CA4
% 
% Created By: Johnathan Tucker
%
% Collaborators: 
%
% The purpose of the script is to act as a driver that will execute the
% functions necesary to solve questions two and three of CA4. Note that the
% solution to question 1 is the PLLT function itself
%
% Created Date: 4/7/2020
%
% Change Log: 
%           - 4/7/2020 Code the PLLT function
%           - 4/8/2020 Code up questions two and three
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
tic
%Global formatting commands to imporve graphing looks:
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex'); 
%% Question 1 Solution
% Please see the attached PLLT equation

%% Question 2 Solution
Question_2();

%% Question 3 Solution
Question_3();

%% Functions Called
% The following functions were built and called as a part of this
% assignment
%
% <include>PLLT.m Vortex_Panel.m Question_2.m Question_3.m NACA_Airfoil.m</include>

toc