%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver script for the fourth question in homework A1
%
% Created by: Johnathan Tucker
%
% Date Created: 10/20/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all;
close all;
clc;
%% Create each rotation matrix
R_1 = [1, 0, 0; 0, cosd(35), sind(35); 0, -sind(35), cosd(35)];
R_2 = [cosd(35), 0, -sind(35); 0, 1, 0; sind(35), 0, cosd(35)];
R_3 = [cosd(35), sind(35), 0; -sind(35), cosd(35), 0; 0, 0, 1];

%% Create each DCM
Q_1 = R_3*R_2*R_1;
Q_2 = R_1*R_2*R_3;
Q_3 = R_3*R_1*R_3;
Q_4 = R_1;
%% Get the eigenstuff for each DCM
[V_1,D_1] = eig(Q_1);
[V_2,D_2] = eig(Q_2);
[V_3,D_3] = eig(Q_3);
[V_4,D_4] = eig(Q_4);
%% Solve for rotation angles
ang_1 = imag(log(D_1(1,1)))*(180/pi);
ang_2 = imag(log(D_2(1,1)))*(180/pi);
ang_3 = imag(log(D_3(2,2)))*(180/pi);
ang_4 = imag(log(D_4(1,1)))*(180/pi);
%% Check rotation angle answer
check_1 = acos(0.5*(trace(Q_1)-1))*(180/pi);
check_2 = acos(0.5*(trace(Q_2)-1))*(180/pi);
check_3 = acos(0.5*(trace(Q_3)-1))*(180/pi);

C_1 = cos(35*pi/180)*eye(3) + (1-cos(35*pi/180))*([1;0;0])*([1;0;0])' - sin(35*pi/180)*[0,0,0;0,0,-1;0,1,0]; 