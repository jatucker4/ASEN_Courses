%% ASEN 3111 - Computational Assignment 3 - Main
%The CA3 code solves questions from three main problem areas: 

%Problem 1:
%Problem 1 requires writing a function to implement the vortex panel 
%method for an arbitrary two dimensional body defined by a set of (x,y) 
%coordinates for boundary points on its surface. The function operates
%effectivly as a translation of the fortran code provided in the Kuethe and
%Chow documentation. 

%Problem 2:
%Problem 2 requires considering a NACA0012 airfoil at alpha = 0 deg using
%varying numbers of vortex panels N to quantitativly assess the error in
%the vortex panel results. A nominal number of panels based on this result
%is then used evaluate the coefficients of lift and pressure for a NACA0012
%at a range of angles of attack. 

%Problem 3: 
%Problem 3 requires assessing plots of c_l vs. alpha for NACA 0012, 2412,
%4412, and 2424 airfoils. These plots are evaluated to consider zero lift
%alpha and lift curve slope approximations in each case, which are then
%compared with thin airfoil theory (TAT). 

% Author: Michael Martinson 
% Collaborators: Fernando Palafox, Jonathan Tucker, Trace Valade
% Date 3/30/2020

%% Housekeeping: 
clc; clear all; close all; 

%Start timer: 
tic

%Global formatting commands:
set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');

%% Problem 1 (See Vortex_Panel.m as solution)

% %Vortex_Panel Function testing operations (Commented out intentionally): 
% %NACA 2412
% m2412 = 2/100;
% p2412 = 4/10; 
% t2412 = 12/100;
% c = 1; %m 
% N = 200; 
% 
% [xtest,ytest] = NACA_Airfoils(m2412,p2412,t2412,c,N);
% 
% %Exact N=12 x y values from the fortran code:
% % xtest = ...
% %     [1 .933 .750 .5 .25 .067 0 .067 .25 .5 .75 .933 1];
% % ytest = ...
% %[0 -0.005 -0.017 -0.033 -0.042 -0.033 0 0.045 0.076 0.072 0.044 0.013 0];
% 
% alpha = 8*(pi/180); %rad
% V_inf = 30; %m/s
% CpSwitch = 0; %If 1, Cp vs. x/c will be plotted:
% 
% %NACA digits for the airfoil considered: 
% Digits2412 = 2412;
% c_l = Vortex_Panel(xtest,ytest,V_inf,alpha,CpSwitch,Digits2412)

%% Problem 2
Nnominal = CA3Prob2Function(); 

%% Problem 3
%Set an N value to use for problem 3 analysis (use 100 fewer panels than 
%the nominal value computed in problem 2 to improve computationanl 
%efficiency with minimal impact on results):
N3 = 500; 

%Call the problem function: 
CA3Prob3Function(N3)

%% end timer: 
fprintf('\n\n'); 
toc

%% Functions Called
% The following functions were built and called as part of this assignment.
%
% <include>Vortex_Panel.m</include>
%
% <include>Prob3NACA_Airfoils.m</include>
%
% <include>NACA_Airfoils.m</include>
%
% <include>CA3Prob3Function.m</include>
%
% <include>CA3Prob2Function.m</include>