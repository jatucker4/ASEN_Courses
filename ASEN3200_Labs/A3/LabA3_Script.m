%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script for lab A3
%
% Created by: Johnathan Tucker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clc;
clear all;
close all;
%% Start script
I = zeros(4,1);
I(1) = dat('Task1T12D2.txt');
I(2) = dat('Task1T8D8.txt');
I(3) = dat('Task1T6D6.txt');
I(4) = dat('Task1T5D5.txt');

T = [12; 8; 6; 5];
D = [2;  8; 6; 5];

T = table(T,D,I);
T.Properties.VariableNames = {'Torque_mNm','Time_s','Moment_of_Inertia_'};
disp(T);

i = 0.00658727;
syms z

zeta = solve(0.1 == exp(-z*pi/sqrt(1-z^2)),z);
zeta = double(zeta);
zeta = zeta(1);
wn   = 3 / (zeta * 1.4065);
wn = 3.6208;
K1 = wn^2;
K2 = 2*zeta*wn;

poles_real = -zeta*wn;
poles_imag = wn*sqrt(zeta^2 - 1);
% Create a matlab simulation for a step response using these gains
n1 = K1;
d1 = 1;
d2 = K2;
d3 = K1;
sysTF = tf(n1,[d1 d2 d3]);
[x,t] = step(sysTF);
x = x*0.5;
unitstep = 0.5*(t>=0);
figure(20)
plot(t,unitstep)
hold on
plot(t,x)
xlabel("$Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
title("$Predicted\:System\:Behavior$",'Interpreter','latex','FontSize',26)
legend("Reference Position", "Actual Position")
%% Pull out the data for question 4
q4_data = readmatrix('k1_49_k2_26_3.csv');
q4_data(1,:) = [];
% q4_data(:,4) = q4_data(:,4).*0;
q4_data(:,1) = q4_data(:,1) - q4_data(1,1) + 1;
q4_data(:,1) = q4_data(:,1) ./ 1000;
q4_data(2689:end,:) = [];
figure(1)
plot(q4_data(:,1),q4_data(:,2))
hold on
plot(q4_data(:,1),q4_data(:,3))
xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
title("$Calculated\:Gains\:K_1\:23.8[mNm/rad]\:and\:K_2\:28.2[mNm/rad/s]$",'Interpreter','latex','FontSize',26)
legend("Reference Position", "Actual Position")

% Plot the torque that the motor creates over time
figure(19)
plot(q4_data(:,1),q4_data(:,4).*33.5)
xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Reaction\:Wheel\:Torque\:[mNm]$",'Interpreter','latex','FontSize',26)
title("$Reaction\:Wheel\:Torque\:Over\:Time$",'Interpreter','latex','FontSize',26)
legend("Step Reference", "Torque Curve")
%% Pull out and process the data for question 5
q5_set_1 = readmatrix('question_5_k1_10.csv');
q5_set_1(1,:) = [];
q5_set_1(:,1) = q5_set_1(:,1) - q5_set_1(1,1) + 1;
q5_set_1(:,1) = q5_set_1(:,1) ./ 1000;

figure(2)
plot(q5_set_1(:,1),q5_set_1(:,2))
hold on
plot(q5_set_1(:,1),q5_set_1(:,3))
xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
title("$Isolated\:Proporional\:Gain\:K_1=10[mNm/rad]\:$",'Interpreter','latex','FontSize',26)
legend("Reference Position", "Actual Position")


q5_set_2 = readmatrix('question_5_k1_15.csv');
q5_set_2(1,:) = [];
q5_set_2(:,1) = q5_set_2(:,1) - q5_set_2(1,1) + 1;
q5_set_2(:,1) = q5_set_2(:,1) ./ 1000;

figure(3)
plot(q5_set_2(:,1),q5_set_2(:,2))
hold on
plot(q5_set_2(:,1),q5_set_2(:,3))
xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
title("$Isolated\:Proporional\:Gain\:K_1=15[mNm/rad]\:$",'Interpreter','latex','FontSize',26)
legend("Reference Position", "Actual Position")



% q5_set_3 = readmatrix('question_5_k1_20.csv');
% q5_set_3(1,:) = [];
% q5_set_3(:,1) = q5_set_3(:,1) - q5_set_3(1,1) + 1;
% q5_set_3(:,1) = q5_set_3(:,1) ./ 1000;
% 
% figure(4)
% plot(q5_set_3(:,1),q5_set_3(:,2))
% hold on
% plot(q5_set_3(:,1),q5_set_3(:,3))
% xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
% ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
% title("$Isolated\:Proporional\:Gain\:K_1=20[mNm/rad]\:$",'Interpreter','latex','FontSize',26)
% legend("Reference Position", "Actual Position")

%% Grab, process, and plot the data for question 6

q6_set_1 = readmatrix('question_6_k2_15.csv');
q6_set_1(1,:) = [];
q6_set_1(:,1) = q6_set_1(:,1) - q6_set_1(1,1) + 1;
q6_set_1(:,1) = q6_set_1(:,1) ./ 1000;

figure(5)
plot(q6_set_1(:,1),q6_set_1(:,2))
hold on
plot(q6_set_1(:,1),q6_set_1(:,3))
xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
title("$Isolated\:Derivative\:Gain\:K_2=15[mNm/rad/s]\:$",'Interpreter','latex','FontSize',26)
legend("Reference Position", "Actual Position")

% q6_set_2 = readmatrix('question_6_k2_30.csv');
% q6_set_2(1,:) = [];
% q6_set_2(9006:end,:) = [];
% q6_set_2(:,1) = q6_set_2(:,1) - q6_set_2(1,1) + 1;
% q6_set_2(:,1) = q6_set_2(:,1) ./ 1000;
% 
% figure(6)
% plot(q6_set_2(:,1),q6_set_2(:,2))
% hold on
% plot(q6_set_2(:,1),q6_set_2(:,3))
% xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
% ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
% title("$Isolated\:Derivative\:Gain\:K_2=30[mNm/rad/s]\:$",'Interpreter','latex','FontSize',26)
% legend("Reference Position", "Actual Position")



q6_set_3 = readmatrix('question_6_k2_35.csv');
q6_set_3(1,:) = [];
q6_set_3(:,1) = q6_set_3(:,1) - q6_set_3(1,1) + 1;
q6_set_3(:,1) = q6_set_3(:,1) ./ 1000;

figure(7)
plot(q6_set_3(:,1),q6_set_3(:,2))
hold on
plot(q6_set_3(:,1),q6_set_3(:,3))
xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
title("$Isolated\:Derivative\:Gain\:K_2=35[mNm/rad/s]\:$",'Interpreter','latex','FontSize',26)
legend("Reference Position", "Actual Position")

%% Grab, process, and plot the data for question 7
q7_data_1 = readmatrix('question_7_k1_49_k2_26_k3_10.csv');
q7_data_1(1,:) = [];
q7_data_1(:,1) = q7_data_1(:,1) - q7_data_1(1,1) + 1;
q7_data_1(:,1) = q7_data_1(:,1) ./ 1000;
q7_data_1(2855:end,:) = [];
figure(8)
plot(q7_data_1(:,1),q7_data_1(:,2))
hold on
plot(q7_data_1(:,1),q7_data_1(:,3))
xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
title("$Calculated\:Gains\:With\:Integral\:Control\:K3=10$",'Interpreter','latex','FontSize',26)
legend("Reference Position", "Actual Position")


q7_data_2 = readmatrix('question_7_k1_49_k2_26_k3_15.csv');
q7_data_2(1,:) = [];
q7_data_2(:,1) = q7_data_2(:,1) - q7_data_2(1,1) + 1;
q7_data_2(:,1) = q7_data_2(:,1) ./ 1000;
q7_data_2(3393:end,:) = [];
figure(9)
plot(q7_data_2(:,1),q7_data_2(:,2))
hold on
plot(q7_data_2(:,1),q7_data_2(:,3))
xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
title("$Calculated\:Gains\:With\:Integral\:Control\:K3=15$",'Interpreter','latex','FontSize',26)
legend("Reference Position", "Actual Position")


q7_data_3 = readmatrix('question_7_k1_49_k2_26_k3_15.csv');
q7_data_3(1,:) = [];
q7_data_3(:,1) = q7_data_3(:,1) - q7_data_3(1,1) + 1;
q7_data_3(:,1) = q7_data_3(:,1) ./ 1000;
q7_data_3(2860:end,:) = [];
figure(10)
plot(q7_data_3(:,1),q7_data_3(:,2))
hold on
plot(q7_data_3(:,1),q7_data_3(:,3))
xlabel("$Arduino\:Time\:[s]$",'Interpreter','latex','FontSize',26)
ylabel("$Integrated\:Gyro\:Output\:[rad]$",'Interpreter','latex','FontSize',26)
title("$Calculated\:Gains\:With\:Integral\:Control\:K3=15$",'Interpreter','latex','FontSize',26)
legend("Reference Position", "Actual Position")


%% Begin function definitions
function inertia = dat(filename)

fileID = fopen(filename);
C = textscan(fileID,'%f %f %f %f');
fclose(fileID);

time    = C{1} * 0.001;     % Convert from ms to s
torque  = C{2} * 0.001;     % Convert from mNm to Nm
w       = C{3} * 2*pi/60;   % Convert from rpm to rad/s

alpha = diff(w) ./ diff(time);  % Derivative of omega wrt time

I = zeros(1,length(alpha));     % Preallocation

for i=1:length(alpha)
    
    if alpha(i)~=0              % Avoid dividing by zero
        I(i) = torque(i+1) / alpha(i);
    end
    
end

inertia = abs(trapz(time(2:end),I));    % abs() b/s inertia>0

end
