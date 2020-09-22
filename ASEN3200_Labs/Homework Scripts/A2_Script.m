%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver script for the fourth question in homework A2
%
% Created by: Johnathan Tucker
%
% Date Created: 10/28/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all;
close all;
clc;
%% Create inertia matrix
I = [783.5, 351.7, 40.27; 351.7, 783.5, -80.27; 40.27, -80.27, 783.5];

%% Calculate eigenstuff
[V, D] = eig(I);

%% Check that the DCM is equal to D
check_mat = ((V')*I*V);
