%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is for problem 21-21 in homework 8
%
% Created by: Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
%% Create Variables
T = 5780;
C_1 = 3.74177e8;
C_2 = 1.43878e4;

% lambda = logspace(-2,-3);
lambda = logspace(-2,3);

for i = 1:length(lambda)
    E(:,i) = C_1/((lambda(i)^5)*(exp(C_2/(lambda(i)*T)) - 1));
end

loglog(lambda,E)
xlabel("$Wavelength\:[\mu m]$",'Interpreter','latex','FontSize',15)
ylabel("$Emissive\:Power\:[W/(m^2 \mu m)]$",'Interpreter','latex','FontSize',15)
ylim([10^-6,10^8])
title("$Wavelength\:vs\:Emissive\:Power$",...
    'Interpreter','latex','FontSize',15)