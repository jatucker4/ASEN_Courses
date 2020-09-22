%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used as the driver to convert the velocity and positions
% to orbital elements, and to envoke Gibbs method to get the velocity and
% position of the NEO.
%
% Created on: 9/17/2019
%
% Created by: Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc; close all; clear all;
%% Begin the first part of the lab
% Create constants
r = [-3424.7; -47.5; 1172];
v = [-0.4250; -3.333; -0.89301];
mu = 42828;
% Calculate orbital elements
orbital_elements_1 = r_v_to_elements(r,v,mu);
% Calculate the orbital period for the propogator
T_1 = sqrt(orbital_elements_1(2)^3/mu)*(2*pi);

%% Second part of the lab
% Create constants
AU = 149597870.7;
r_1 = [0.5887,-0.2206,0.0239] * AU;
r_2 = [0.5027,0.2289,0.0436] * AU;
r_3 = [0.3243,0.4560,0.0453] * AU;
mu_sun = 132712000000;
% Get the velocity corresponding to the second position vector using Gibb's
% method.
[r,v] = Gibbs_method(r_1,r_2,r_3,mu_sun);
% Get the orbital elements from the position and velocity vectors
orbital_elements_2 = r_v_to_elements(r,v,mu_sun);
% Get the period for the propogator.
T_2 = sqrt(orbital_elements_2(2)^3/mu_sun)*(2*pi);

r_p = orbital_elements_2(2)*(1-orbital_elements_2(3));
r_a = 2*orbital_elements_2(2) - r_p;