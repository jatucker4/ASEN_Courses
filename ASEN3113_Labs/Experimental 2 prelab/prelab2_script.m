%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for the second pre-lab
% 
% Created by: Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
close all;
clear all;
%% Solve for slope and inital temp
% Create vectors
d_1 = 1 + 3/8;
d_vec = [0, d_1, d_1+0.5 , d_1+1, d_1+1.5, d_1+2, d_1+2.5, d_1+3, d_1+3.5];
T_vec = [18.53, 22.47, 26.87, 30.05, 35.87, 38.56, 41.50, 46.26];
% Apply polyfit
p = polyfit(d_vec(2:end),T_vec,1);
% This means the line looks like y = p(1)x + p(2) so to get y(0)
T_0 = p(2);
H = p(1);
L = d_vec(9) + 1;
%% Create heat eqn at the last thermocouple
temp = -8*7.8607*5.875;
temp2 = (pi*4.875)/(2*5.875);   
temp3 = (pi)/(2*5.875);
t = 1;
x_var = 0:10;
syms n
increment = 1:10;
u_2 = zeros(1,length(n)+1);
u_2(1) = (7.949 + 7.8607*4.875);
for i = 1:length(increment)
    u_temp = (7.949 + 7.8607*4.875) + symsum(((temp*((-1)^(n + 1)))/(((2*n - 1)*pi)^2))* ...
        sin((2*n - 1)*temp2) * exp(-(((2*n - 1)*temp3)^2) * .07409 * t) , n, 1, increment(i));
    u_2(i+1) = double(u_temp);
end
figure(1)
plot(x_var,u_2,'o')
xlabel('Interval Value','Interpreter','latex','FontSize',26)
ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
ylim([0 50])
title("Temperature at t=1 Versus Index",'Interpreter','latex','FontSize',26)

t = 1000;
syms n
increment = 1:10;
u_2 = zeros(1,length(n)+1);
u_2(1) = (7.949 + 7.8607*4.875);
for i = 1:length(increment)
    u_temp = (7.949 + 7.8607*4.875) + symsum(((temp*((-1)^(n + 1)))/(((2*n - 1)*pi)^2))* ...
        sin((2*n - 1)*temp2) * exp(-(((2*n - 1)*temp3)^2) * .07409 * t) , n, 1, increment(i));
    u_2(i+1) = double(u_temp);
end
figure(2)
plot(x_var,u_2,'o')
xlabel('Interval Value','Interpreter','latex','FontSize',26)
ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
ylim([0 50])
title("Temperature at t=1000 Versus Index",'Interpreter','latex','FontSize',26)
%% Question 6
t = 1:1000;
syms n
increment = 1:10;
u_2 = zeros(1,length(n)+1);
% u_2(1) = (7.949 + 7.8607*4.875);
for i = 1:length(t)
    u_2(i) = (7.949 + 7.8607*4.875) + symsum(((temp*((-1)^(n + 1)))/(((2*n - 1)*pi)^2))* ...
        sin((2*n - 1)*temp2) * exp(-(((2*n - 1)*temp3)^2) * .07409 * t(i)) , n, 1, 1);
%     u_2(i+1) = double(u_temp);
end
figure(3)
plot(t,u_2)
xlabel('Time[s]','Interpreter','latex','FontSize',26)
ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
ylim([0 50])
title("Temperature at Last Thermocouple Versus Time",'Interpreter','latex','FontSize',26)

t = 1:1000;
alpha_vec = [0.25, 0.5, 1, 1.25, 1.5]*.07409;
syms n
for j = 1:length(alpha_vec)
for i = 1:length(t)
    u_2(i) = (7.949 + 7.8607*4.875) + symsum(((temp*((-1)^(n + 1)))/(((2*n - 1)*pi)^2))* ...
        sin((2*n - 1)*temp2) * exp(-(((2*n - 1)*temp3)^2) * alpha_vec(j) * t(i)) , n, 1, 1);
end
figure(4)
plot(t,u_2)
hold on
xlabel('Time[s]','Interpreter','latex','FontSize',26)
ylabel("$Temperature\:[C]$",'Interpreter','latex','FontSize',26)
ylim([0 50])
title("Temperature at Last Thermocouple Versus Time",'Interpreter','latex','FontSize',26)
end
legend("0.25\alpha","0.5\alpha","\alpha","1.25\alpha","1.5\alpha")