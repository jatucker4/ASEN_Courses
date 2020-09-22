function [Me_out,Pe] = guess_iteration(A2oAt)
%GUESS_ITERATION Summary of this function goes here
%   Detailed explanation goes here
% Get M1
M1 = eq10_23_solve_M(A2oAt,1);
% Get M2
M2 = sqrt((1 + (.4/2)*M1^2)/(1.4*M1^2 - .4/2));
% Get change in enthalpy
enth = abs(get_delta_s(M2));
% Get p02op01
p02op01 = exp(-enth/287);
% Get A2/Astar
A2oAstar = eq10_23_solve_area(M2);
% Get Ae/Astar
AeoAstar_c2 = 1.53*(A2oAt^-1)*A2oAstar;
% Solve for Me
Me_out = eq10_23_solve_M(AeoAstar_c2,0);
% Solve Po2ope
po2ope = (1+(.4/2)*Me_out^2)^(1.4/.4);
Pe = (po2ope^-1)*(p02op01);
end

