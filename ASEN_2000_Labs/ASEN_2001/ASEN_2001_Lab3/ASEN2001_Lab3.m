%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2001 Lab3
% Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear;
close all;
clc;
%% Variables
foam_modulus = 0.035483*(10^9);
balsa_modulus = 3.2953*(10^9);
[NUM, ~, ~] = xlsread('TestData.xlsx');
test_number = NUM(:,1);
force_failure = NUM(:,2);
a = NUM(:,3);
w = NUM(:,4);
d_f = NUM(:,5);
t_balsa = 0.03125/39.37;
t_foam = 0.75/39.37;
t_total = 2*t_balsa + t_foam;
L = 36/39.37;
b = L - 2*a;
%% Calculating Reaction forces for each test
reaction_b = (((force_failure*0.5).*a) + (force_failure*0.5).*(a+b))/L;
reaction_a = (force_failure*0.5) + (force_failure*0.5) - reaction_b;
%% Segment One shear and moments
% Create x-matrix
for i = 1:length(a)
    temp = a(i,:);
    x(i,:) = linspace(0,temp,24);
end
shear_segment_one = reaction_a;
for i = 1:length(a)
    moment_segment_one(i,:) = shear_segment_one'.*x(i,:);
end
%% Segment Two shear and moments
% Create x_two matrix
for i = 1:24
    x_two(i,:) = linspace(0,(L/2),24);
end
shear_segment_two = reaction_a - force_failure;
for i = 1:24
    moment_segment_two(i,:) = shear_segment_two'.*x_two(i,:) + (force_failure'/2).*a(i,:)';
end

%% Combine the moment matrices
moment_total = [moment_segment_one , moment_segment_two];
max_moments = max(moment_total,[],2);
%% Calculate the stress at failure for each test
[I_b, I_f] = beammoment(w);
c = t_balsa + (0.5*t_foam);
sigma_fail = (-max_moments*c)./(I_b + (foam_modulus/balsa_modulus)*I_f);
%% Now the shear of failue
shear_total = [shear_segment_one, shear_segment_two];
max_shear = max(shear_total, [], 2);
A_f = w*t_total;
shear_fail = (3/2)*(max_shear./A_f);
%% Eliminating Outliers
outliers = [2,7,11,16,18,19, 24];
sigma_fail(outliers) = [];
sigma_fail = mean(sigma_fail);
shear_fail(outliers) = [];
shear_fail = mean(shear_fail);
%% Shear and moment diagrams
length_1 = 3*a(1,1);
x_coords_moment = linspace(0,length_1,48);
moment_plot_1 = sort(moment_segment_one(1,:),'ascend');
moment_plot_2 = sort(moment_segment_two(1,:),'descend');
moment_plot = [moment_plot_1,moment_plot_2];
figure(1)
p(1) = plot(x_coords_moment,moment_plot,'DisplayName','Moment');
title('$Moment\:Versus\:X\:Position$','Interpreter','latex')
xlabel('$X\:(m)$','Interpreter','latex')
ylabel('$Moment\:(Nm)$','Interpreter','latex')
legend([p(1)])
shear_plot = [shear_segment_one(1,1),shear_segment_two(1,1)];
x_coords_moment = [0,length_1];
figure(2)
p(2) = plot(x_coords_moment,shear_plot,'DisplayName','Shear');
title('$Shear\:Versus\:X\:Position$','Interpreter','latex')
xlabel('$X\:(m)$','Interpreter','latex')
ylabel('$Shear\:(N)$','Interpreter','latex')
legend([p(2)])