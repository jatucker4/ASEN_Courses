%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3112 - Homework 7
% 
% Created By: Johnathan Tucker
%
% Collaborators: N/A
%
% The purpose of the script is to perform calculations for the DSM-FEM
% problem(7.2) on homework 7
%
% Created Date: 3/14/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
close all;
clear all;
%% Create necessary variables
E = 3000;
L_1 = 30;
L_2 = 20*sqrt(3);
L_3 = 20;
A_1 = 2;
A_2 = 4;
A_3 = 3;
K_hat_1 = [200, 0 ; 0 , 0];
%% Solve for Global K matrix for bar 2
phi_2 = -30;
K_hat_2 = (E*A_2/L_2).*[(cosd(phi_2))^2 , sind(phi_2)*cosd(phi_2) ; sind(phi_2)*cosd(phi_2) , (sind(phi_2)^2)];
K_mat_2 = [K_hat_2 , -K_hat_2 ; -K_hat_2, K_hat_2];
fprintf("Global K2 matrix is:\n")
disp(K_mat_2)
fprintf("\n")
%% Solve for Global K matrix for bar 3
phi_3 = -120;
K_hat_3 = (E*A_3/L_3).*[(cosd(phi_3))^2 , sind(phi_3)*cosd(phi_3) ; sind(phi_3)*cosd(phi_3) , (sind(phi_3)^2)];
K_mat_3 = [K_hat_3 , -K_hat_3 ; -K_hat_3, K_hat_3];
fprintf("Global K3 matrix is:\n")
disp(K_mat_3)
fprintf("\n")
%% Merge the global K matrices
K_global = [K_hat_1 , zeros(2,2) , zeros(2,2) , -K_hat_1 ; ...
    zeros(2,2) , K_hat_2 , zeros(2,2) , -K_hat_2;...
    zeros(2,2), zeros(2,2), K_hat_3, -K_hat_3;...
    -K_hat_1, -K_hat_2, -K_hat_3, K_hat_1+K_hat_2+K_hat_3];
fprintf("Full global K matrix is:\n")
disp(K_global)
fprintf("\n")
%% Apply Boundary Conditions
f = [0; -200];
K_global_reduced = K_hat_1+K_hat_2+K_hat_3;
fprintf("Reduced K matrix is:\n")
disp(K_global_reduced)
fprintf("\n")
%% Solve for the displacements
u_vec = K_global_reduced\f;
fprintf("Displacement solutions are ux4 = %f and uy4 = %f\n",u_vec(1,1), u_vec(2,1))
%% Solve for reaction forces
u_vec_global = [zeros(6,1) ; u_vec];
f_reactions = K_global*u_vec_global;
fprintf("Reaction forces are:\n")
disp(f_reactions)
fprintf("\n")
%% Solve for internal forces in the bars
F_1 = f_reactions(1,1)*cos(0);
F_2 = f_reactions(3,1)*cosd(phi_2) + f_reactions(4,1)*sind(phi_2);
F_3 = f_reactions(5,1)*cosd(phi_3) + f_reactions(6,1)*sind(phi_3);
fprintf("Internal forces in the bars are: F1 = %f, F2 = %f, F3 = %f\n",F_1,F_2,F_3)
