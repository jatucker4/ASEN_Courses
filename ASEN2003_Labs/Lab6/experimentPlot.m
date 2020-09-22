%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003: Lab 6 Rotary Position Control
%   Group 13:
%   Sam Hartman
%   Joshua Seedorf
%   Johnathan Tucker
%   Sean Yoo
%   
% The purpose of this script is to plot the behavior of both the
% rigid and the flexible arm from the experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

%% Read in data
raw = load('modeled_gains_K1_50_K2_1.5');
%Identify the start of a step
theta_d = 0.379;
% Pull necessary data
start = find(raw(:,6)>abs(theta_d*0.9));
data = raw(start(1):size(raw,1),:);
t       = (data(:,1) - data(1,1))/1000;   %[s]
theta   = data(:,2);        %[rad]
d       = data(:,3);        %[m]
alpha   = data(:,4);        %[rad/s]
V_tip   = data(:,5);        %[m/s]
voltage = data(:,7);        %[V]
K1      = data(1,8);
K2      = data(1,9);
K3      = data(1,10);
K4      = data(1,11);
if K2 == 0 && K4 ==0
    [zeta,w_n] = getZeta(K1,K3);
    t_s = -log(0.05)./(zeta*w_n);
else
    %% Find settling time
    cross = find(theta > 0.95*theta_d,1, 'last');
    t_s = t(cross(1));
end

%% Plot
theta = smooth(theta);
%Theta VS Theta_d
figure()
hold on
plot([0 max(t)],[theta_d theta_d])
plot(t,theta)
plot([0 max(t)],1.05*theta_d*[1 1],'k--')
plot([0 max(t)],0.95*theta_d*[1 1],'k--')

legend('Reference','Actual','5% Margin')
ylabel('$\theta\:[rad]$', 'interpreter','latex','FontSize',20)
xlabel('$Time\:[s]$', 'interpreter','latex','FontSize',20)
title('$Arm\:Angular\:Position\:Over\:Time\:(Experimental)$',...
    'Interpreter','latex','FontSize',20)
xlim([0 1])

%Voltage
figure()
hold on
plot(t,voltage)
ylabel('Voltage [V]')
xlabel('time [s]')
title('Voltage Response')
xlim([0 5])

% Deflection
figure()
hold on
% plot([0 max(t)],[theta_d theta_d])
d = smooth(d,'rloess');
plot(t,d)
legend('Tip displacement')
ylabel('$Tip\:Deflection\:[m]$', 'interpreter','latex','FontSize',20)
xlabel('$Time\:[s]$', 'interpreter','latex','FontSize',20)
title('$Arm\:Tip\:Deflection\:Over\:Time\:(Experimental)$',...
    'Interpreter','latex','FontSize',20)
xlim([0 1])