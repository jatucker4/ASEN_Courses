%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver script for the fourth question in homework A6
%
% Created by: Johnathan Tucker
%
% Date Created: 12/9/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all;
close all;
clc;
%% Create the rotation matrix
Q_a = [0, 0, 1; -1, 0, 0; 0, -1, 0];
Q_b = [0, 1, 0; -1/sqrt(2), 0, 1/sqrt(2); 1/sqrt(2), 0, 1/sqrt(2)];
Q_c = [-1, 0, 0; 0, -1/sqrt(2), 1/sqrt(2); 0, 1/sqrt(2), 1/sqrt(2)];
[V_1,D_1] = eig(Q_a);
[V_2,D_2] = eig(Q_b);
[V_3,D_3] = eig(Q_c);
%% Solve for rotation angles
ang_1 = imag(log(D_1(1,1)))*(180/pi);
ang_2 = imag(log(D_2(1,1)))*(180/pi);
ang_3 = imag(log(D_3(2,2)))*(180/pi);
%% Check rotation angle answer
check_1 = acos(0.5*(trace(Q_a)-1))*(180/pi);
check_2 = acos(0.5*(trace(Q_b)-1))*(180/pi);
check_3 = acos(0.5*(trace(Q_c)-1))*(180/pi);