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
R_1 = [1, 0, 0; 0, cosd(60), sind(60); 0, -sind(60), cosd(60)];
R_2 = [cosd(60), 0, -sind(60); 0, 1, 0; sind(60), 0, cosd(60)];
R_3 = [cosd(60), sind(60), 0; -sind(60), cosd(60), 0; 0, 0, 1];

%% Create each DCM
Q_1 = R_3*R_2*R_1;
Q_2 = R_1*R_2*R_3;
Q_3 = R_3*R_1*R_3;

%% Get the eigenstuff for each DCM
[V_1,D_1] = eig(Q_1);
[V_2,D_2] = eig(Q_2);
[V_3,D_3] = eig(Q_3);

%% Solve for rotation angles
ang_1 = imag(log(D_1(2,2)))*(180/pi);
ang_2 = imag(log(D_2(1,1)))*(180/pi);
ang_3 = imag(log(D_3(2,2)))*(180/pi);

%% Check rotation angle answer
check_1 = acos(0.5*(trace(Q_1)-1))*(180/pi);
check_2 = acos(0.5*(trace(Q_2)-1))*(180/pi);
check_3 = acos(0.5*(trace(Q_3)-1))*(180/pi);

%% Get unit quaternion for each
q_1 = [sind(check_1/2)*V_1(:,1); cosd(check_1/2)];
q_2 = [sind(check_2/2)*V_2(:,3); cosd(check_2/2)];
q_3 = [sind(check_3/2)*V_3(:,1); cosd(check_3/2)];