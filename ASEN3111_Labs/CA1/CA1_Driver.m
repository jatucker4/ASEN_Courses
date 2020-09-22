%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3111 - CA1 
% 
% Created By: Johnathan Tucker
%
% Collaborators: N/A
%
% The purpose of the script is to contain all of the constants and basic
% computation that will be passed to the question functions where all of
% the necessary outputs for each question will come from
%
% Created Date: 1/21/2020
%
% Change Log: 
%           - 1/25/2020: Code up trap and simpsons rule for the cylinder
%           - 1/28/2020: Code up trap rule and error logic for the airfoil
%           - 2/02/2020: Finalize formatting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
tic
%% Create all necessary constants and put them in a struct
const.rho_inf = 1.225; % kg/m^3
const.p_inf = 101.3e3; % Pa
const.d = 1; % m
const.r = const.d/2; % m
const.v_inf = 30; % m/s
const.q_inf = 0.5*const.rho_inf*const.v_inf^2; % Pa
const.gamma = 2*pi*const.r*const.v_inf; % m^2/s^2
const.aoa = 9; % Degrees
const.chord = 2; % m

%% Begin Question 1
Question_1(const);
%% Begin Question 2
Question_2(const);
toc