%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the driver script for the fourth question in homework 3
%
% Created by: Johnathan Tucker
%
% Date Created: 9/17/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all;
close all;
clc;
%% Create all needed variables here
r = [4973; -1798; -1748];
v = [-0.5713; 0.4174; 2.136];
k = [0, 0, 1];
mu = 22030;
% Use above to get ang momentum, eccentricity and n
h = cross(r,v);
e = (cross(v,h)/mu) - (r/norm(r));
n = cross(k,h);
% Get the norms of each
e_norm = norm(e);
h_norm = norm(h);
n_norm = norm(n);
r_norm = norm(r);
v_norm = norm(v);

a = ((h_norm^2)/mu)/(1-e_norm^2);
% Calculate the angles with the checks
i = acos(h(3)/h_norm);
if(n(2) <0)
    cap_ohm = (2*pi) - acos(n(1)/n_norm);
else
    cap_ohm = acos(n(1)/n_norm);
end

if(e(3) <0)
    ohm = (2*pi) - acos(dot(n,e)/(n_norm*e_norm));
else
    ohm = acos(dot(n,e)/(n_norm*e_norm));
end

if(dot(r,v) < 0)
    theta = (2*pi) - acos(dot(e,r)/(r_norm*e_norm));
else
    theta = acos(dot(e,r)/(r_norm*e_norm));
end
%% Create the DCM

dcm_1 = [cos(ohm), sin(ohm), 0 ; -sin(ohm), cos(ohm), 0; 0, 0, 1 ];
dcm_2 = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
dcm_3 = [cos(cap_ohm), sin(cap_ohm), 0 ; -sin(cap_ohm), cos(cap_ohm), 0; 0, 0, 1 ];


step_1 = dcm_3*v;
step_2 = dcm_2*step_1;
v_peri_2 = dcm_1*step_2;

dcm_final = dcm_1*dcm_2*dcm_3;
v_peri_3 = dcm_final*v;
v_peri_norm = norm(v_peri_2);

fpa = atan((e_norm*sin(theta))/(1+e_norm*cos(theta)));
%% Calculate the perifocal frame position and velocity
r_peri = dcm_final*r;
r_peri_norm = norm(r_peri);
T = sqrt(a^3/mu)*(2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical solution
v_p = (mu/h_norm)*(-sin(theta));
v_q = (mu/h_norm)*(e_norm+cos(theta));