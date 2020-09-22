%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003 Lab 1
%
% Created by: Johnathan Tucker
% This script is the driver script that calls all necessary functions to  
% build the plot of the roller coaster as well as the G plot and speed plot
%
% Last Updated: 1/28/2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear;
close all;
clc;
%% Plot the roller coaster getting the length for each component

% Initial setup constants for the helix
t_start = 0;
r = 25;
loops = 3;
helix_height = 30;
position = [0,0,125];
tot_length = 0;

% Create the figure for the roller coaster path to pass to functions
coaster_plot = figure('Name', 'Path of Roller Coaster');

% Create the helix and plot it
[pos_final, helix_length, helix_time, helix_z_t] = Lab1_Helix(r,...
    loops, helix_height, t_start, position, coaster_plot);
% Get the helix length and append it
tot_length = tot_length + helix_length;
helix = tot_length;
hold on
scatter3(pos_final(1), pos_final(2), pos_final(3), 'filled')
hold on

% Create the first transition
[pos_final_hill, hill_length, hill_time,hill_1_z_t,...
    hill_1_k] = transition(pos_final,coaster_plot,0);
% Get the transitions length
tot_length = tot_length + hill_length;
hill_1 = tot_length;
hold on
scatter3(pos_final_hill(1), pos_final_hill(2), pos_final_hill(3), 'filled')
hold on

% Create the second transition
[pos_final_hill_2, hill_length_2, hill_time_2,hill_2_z_t,...
    hill_2_k] = transition(pos_final_hill, coaster_plot,1);
hold on
scatter3(pos_final_hill_2(1), pos_final_hill_2(2),...
    pos_final_hill_2(3), 'filled')
% Get the length and append it
tot_length = tot_length + hill_length_2;
hill_2 = tot_length;
hold on

% Create the loop
[pos_final_loop, loop_length, loop_time,loop_z_t,...
    loop_gs] = Lab1_Loop(20, 0, pos_final_hill_2, coaster_plot);
% Get the length and append it
tot_length = tot_length + loop_length;
loop = tot_length;
hold on

% Create the bank
[pos_bank_final, bank_length, bank_time,bank_z_t,...
    bank_gs] = Lab1_Bank(36,pos_final_loop,coaster_plot);
hold on
scatter3(pos_bank_final(1),pos_bank_final(2),pos_bank_final(3),'filled')
% Get the length and append it
tot_length = tot_length + bank_length;
bank = tot_length;
hold on

% Create the top hat segment
[parabola_length, pos_final_parabola, parabola_time, parabola_z_t,...
    parabola_gs] = Lab1_Hill(coaster_plot, pos_bank_final);
% Get the length and append it
tot_length = tot_length + parabola_length;
parabola = tot_length;
hold on

% Create the braking segment
[brake_length, brake_pos_final,...
    brake_z_t] = Lab1_Brake(pos_final_parabola,coaster_plot);
% Get the length and append it
tot_length = tot_length + brake_length;
brake = tot_length;
% Display total length
fprintf('Total track length is: %f\n', tot_length)

%% Plot the speed as a function of position
% Create Speed functions
syms z s
speed = sqrt(2*9.81*(125-z));
% Speed function for braking is separate to account for friction
speed_2 = sqrt(49.522722^2 + 2*(-16.35)*(s - parabola));

% Create the length vectors
helix_vec = linspace(0, helix);
hill_1_vec = linspace(helix,hill_1);
hill_2_vec = linspace(hill_1,hill_2);
loop_vec = linspace(hill_2,loop);
bank_vec = linspace(loop,bank);
parabola_vec = linspace(bank,parabola,length(parabola_z_t));
brake_vec = linspace(parabola,brake);

% Create the segment speed vectors
helix_speed = double(subs(speed,z,helix_z_t));
hill_1_speed = double(subs(speed,z,hill_1_z_t));
hill_2_speed = double(subs(speed,z,hill_2_z_t));
loop_speed = double(subs(speed,z,loop_z_t));
bank_speed = double(subs(speed,z,bank_z_t));
parabola_speed = double(subs(speed,z,parabola_z_t));
brake_speed = double(subs(speed_2,s,brake_vec));

% Create the speed plot
figure(2)
plot(helix_vec,helix_speed)
hold on
plot(hill_1_vec,hill_1_speed)
hold on
plot(hill_2_vec,hill_2_speed)
hold on
plot(loop_vec,loop_speed)
hold on
plot(bank_vec,bank_speed)
hold on
plot(parabola_vec,parabola_speed)
hold on
plot(brake_vec,brake_speed)
title('$Speed\:of\:Rollercoaster\:vs\:Position$','Interpreter','latex',...
    'FontSize',16)
xlabel('$Position\:(m)$','Interpreter','latex','FontSize',16)
ylabel('$Speed\:(m/s)$','Interpreter','latex','FontSize',16)
legend('Helix Speed', 'Transition 1 Speed', 'Transition 2 Speed',...
    'Loop Speed', 'Bank Speed', 'Tower Speed', 'Braking speed')
%% G's as a function of S
% Calculate the G's for the helix, transitions, and braking segments
helix_angle = atan(pi*2*r/helix_height);
helix_normal = ones(1,100)*sin(helix_angle);
helix_lateral = helix_speed.^2 / (9.81 * r);
ramp_angle = acos(dot(pos_final_hill_2,pos_final_hill)./...
    (norm(pos_final_hill_2)*norm(pos_final_hill)));
hill_1_gs = ones(1,100)*sin(ramp_angle);
hill_2_gs = hill_2_speed.^2*(hill_2_k/9.81) + 1;
brake_gs = ones(1,100)*(16.35/9.81);
bank_gs = ones(1,100)*bank_gs;

% Create G limits
up_limit = linspace(6,6);
forward_limit = linspace(5,5);
back_limit = linspace(4,4);
down_limit = linspace(1,1);
lateral_limit = linspace(3,3);
x_axis = linspace(0,tot_length);
max_lateral = max([helix_lateral,bank_gs]);
max_normal = max([helix_normal,hill_1_gs,hill_2_gs,loop_gs,parabola_gs]);
max_back = max(brake_gs);
% Create G plot
figure(3)
plot(helix_vec,helix_normal)
hold on
plot(helix_vec,helix_lateral)
hold on
plot(hill_1_vec,hill_1_gs)
hold on
plot(hill_2_vec,hill_2_gs)
hold on
plot(loop_vec,loop_gs)
hold on
plot(bank_vec,bank_gs)
hold on
plot(parabola_vec,parabola_gs)
hold on
plot(brake_vec,brake_gs)
hold on
plot(x_axis,up_limit)
hold on
plot(x_axis,back_limit)
hold on
plot(x_axis,forward_limit)
hold on
plot(x_axis,down_limit)
hold on
plot(x_axis,lateral_limit)
hold on
vline(helix)
hold on
vline(hill_1)
hold on
vline(hill_2)
hold on
vline(loop)
hold on
vline(bank)
hold on
vline(parabola)
hold on
vline(brake)
title('$Gs\:of\:Rollercoaster\:vs\:Position$','Interpreter','latex',...
    'FontSize',16)
xlabel('$Position\:(m)$','Interpreter','latex','FontSize',16)
ylabel('$Gs\:(a/g)$','Interpreter','latex','FontSize',16)
legend('Helix Upward Gs','Helix Lateral Gs', 'Transition 1 Upward Gs',...
    'Transition 2 Upward Gs','Loop Upward Gs',...
    'Banked Turn Lateral Gs','Tower Upward Gs', 'Braking Backward Gs',...
    'Upward Gs Limit', 'Backward Gs Limit', 'Forward Gs Limit',...
    'Downward Gs Limit','Lateral Gs Limit')
%% Guesstimating the time
% Get the average speed
speed_vec = [helix_speed,hill_1_speed,hill_2_speed,...
    loop_speed,parabola_speed,bank_speed,brake_speed];
speed_mean = mean(speed_vec);
% Find time estimate
time_spent = tot_length/speed_mean;
fprintf('The roller coaster takes about %f seconds\n',time_spent)