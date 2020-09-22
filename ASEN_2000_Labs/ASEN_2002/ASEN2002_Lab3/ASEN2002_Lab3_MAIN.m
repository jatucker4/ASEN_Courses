clear all;
clear fig;
clc;
%% position stuff
aoa = -15:16;
port_number = linspace(1,20,20);
x = [0, 0.175, 0.35, 0.7, 1.05, 1.4, 1.75, 2.1,...
                       2.8, 3.5,2.8, 2.1,...
                       1.4, 1.05, 0.7, 0.35, 0.175];
y = [0.14665, 0.33075, 0.4018, 0.476, 0.49, 0.4774, ...
                       0.4403, 0.38325, 0.21875,0, 0, 0, ...
                       0, 0, 0.0014, 0.0175, 0.03885];
%% Retrieving Cp data for every file
cpData_1 = getcpData('AirfoilPressure_S011_G01.csv');
cpData_2 = getcpData('AirfoilPressure_S012_G03.csv');
cpData_3 = getcpData('AirfoilPressure_S012_G05.csv');
cpData_4 = getcpData('AirfoilPressure_S012_G07.csv');
cpData_5 = getcpData('AirfoilPressure_S012_G09.csv');
cpData_6 = getcpData('AirfoilPressure_S012_G11.csv');
cpData_7 = getcpData('AirfoilPressure_S012_G13.csv');
cpData_8 = getcpData('AirfoilPressure_S012_G15.csv');
%% Getting C_l data and organizing
C_l_1_true = [];
C_l_2_true = [];
C_l_3_true = [];
% group 1
[c_l_1,c_l_2,c_l_3,c_d_1,c_d_2,c_d_3] = getC_l(cpData_1,x,y);
% group 3
[c_l_4,c_l_5,c_l_6,c_d_4,c_d_5,c_d_6] = getC_l(cpData_2,x,y);
% group 5
[c_l_7,c_l_8,c_l_9,c_d_7,c_d_8,c_d_9] = getC_l(cpData_3,x,y);
% group 7
[c_l_10,c_l_11,c_l_12,c_d_10,c_d_11,c_d_12] = getC_l(cpData_4,x,y);
% group 9
[c_l_13,c_l_14,c_l_15,c_d_13,c_d_14,c_d_15] = getC_l(cpData_5,x,y);
% group 11
[c_l_16,c_l_17,c_l_18,c_d_16,c_d_17,c_d_18] = getC_l(cpData_6,x,y);
% group 13
[c_l_19,c_l_20,c_l_21,c_d_19,c_d_20,c_d_21] = getC_l(cpData_7,x,y);
% group 15
[c_l_22,c_l_23,c_l_24,c_d_22,c_d_23,c_d_24] = getC_l(cpData_8,x,y);

%% Plots
C_l_1_true_2 = [c_l_22(1,1),c_l_19(1,1),c_l_16(1,1),c_l_13(1,1),c_l_10(1,1),...
    c_l_7(1,1),c_l_4(1,1),c_l_1(1,1),c_l_22(1,2),c_l_19(1,2),c_l_16(1,2),...
    c_l_13(1,2),c_l_10(1,2), c_l_7(1,2),c_l_4(1,2),c_l_1(1,2),...
    c_l_22(1,3),c_l_19(1,3),c_l_16(1,3),c_l_13(1,3),c_l_10(1,3),...
    c_l_7(1,3),c_l_4(1,3),c_l_1(1,3),c_l_22(1,4),c_l_19(1,4),c_l_16(1,4),...
    c_l_13(1,4),c_l_10(1,4), c_l_7(1,4),c_l_4(1,4),c_l_1(1,4)];


C_l_2_true_2 = [c_l_23(1,1),c_l_20(1,1),c_l_17(1,1),c_l_14(1,1),c_l_11(1,1),...
    c_l_8(1,1),c_l_5(1,1),c_l_2(1,1),c_l_23(1,2),c_l_20(1,2),c_l_17(1,2),...
    c_l_14(1,2),c_l_11(1,2), c_l_8(1,2),c_l_5(1,2),c_l_2(1,2),...
    c_l_23(1,3),c_l_20(1,3),c_l_17(1,3),c_l_14(1,3),c_l_11(1,3),...
    c_l_8(1,3),c_l_5(1,3),c_l_2(1,3),c_l_23(1,4),c_l_20(1,4),c_l_17(1,4),...
    c_l_14(1,4),c_l_11(1,4), c_l_8(1,4),c_l_5(1,4),c_l_2(1,4)];

C_l_3_true_2 = [c_l_24(1,1),c_l_21(1,1),c_l_18(1,1),c_l_15(1,1),c_l_12(1,1),...
    c_l_9(1,1),c_l_6(1,1),c_l_3(1,1),c_l_24(1,2),c_l_21(1,2),c_l_18(1,2),...
    c_l_15(1,2),c_l_12(1,2), c_l_9(1,2),c_l_6(1,2),c_l_3(1,2),...
    c_l_24(1,3),c_l_21(1,3),c_l_18(1,3),c_l_15(1,3),c_l_12(1,3),...
    c_l_9(1,3),c_l_6(1,3),c_l_3(1,3),c_l_24(1,4),c_l_21(1,4),c_l_18(1,4),...
    c_l_15(1,4),c_l_12(1,4), c_l_9(1,4),c_l_6(1,4),c_l_3(1,4)];

% NACA data
naca_cl = [-.13,.18,.31,.49,.62,.78,1.08,1.32];
naca_aoa = [-8,-4,-2,0,2,4,8,12];
figure(1)
p(1) = plot(aoa,C_l_1_true_2(1,:),':ro','DisplayName','C_l at 9 m/s');
hold on
p(2) = plot(aoa,C_l_2_true_2(1,:),':bo','DisplayName','C_l at 17 m/s');
hold on 
p(3) = plot(aoa,C_l_3_true_2(1,:),':go','DisplayName','C_l at 34 m/s');
hold on
p(4) = plot(naca_aoa,naca_cl,':ko','DisplayName','NACA ClarkY-14 Data');
title('$Coefficient\:of\:Lift\:Versus\:Angle\:of\:Attack\:$',...
    'Interpreter','latex')
xlabel('$Angle\:of\:Attack(degrees)$','Interpreter','latex')
ylabel('$C_l$','Interpreter','latex')
legend([p(1),p(2),p(3),p(4)])
%% Now C_d
C_d_1_true_2 = [c_d_22(1,1),c_d_19(1,1),c_d_16(1,1),c_d_13(1,1),c_d_10(1,1),...
    c_d_7(1,1),c_d_4(1,1),c_d_1(1,1),c_d_22(1,2),c_d_19(1,2),c_d_16(1,2),...
    c_d_13(1,2),c_d_10(1,2), c_d_7(1,2),c_d_4(1,2),c_d_1(1,2),...
    c_d_22(1,3),c_d_19(1,3),c_d_16(1,3),c_d_13(1,3),c_d_10(1,3),...
    c_d_7(1,3),c_d_4(1,3),c_d_1(1,3),c_d_22(1,4),c_d_19(1,4),c_d_16(1,4),...
    c_d_13(1,4),c_d_10(1,4), c_d_7(1,4),c_d_4(1,4),c_d_1(1,4)];


C_d_2_true_2 = [c_d_23(1,1),c_d_20(1,1),c_d_17(1,1),c_d_14(1,1),c_d_11(1,1),...
    c_d_8(1,1),c_d_5(1,1),c_d_2(1,1),c_d_23(1,2),c_d_20(1,2),c_d_17(1,2),...
    c_d_14(1,2),c_d_11(1,2), c_d_8(1,2),c_d_5(1,2),c_d_2(1,2),...
    c_d_23(1,3),c_d_20(1,3),c_d_17(1,3),c_d_14(1,3),c_d_11(1,3),...
    c_d_8(1,3),c_d_5(1,3),c_d_2(1,3),c_d_23(1,4),c_d_20(1,4),c_d_17(1,4),...
    c_d_14(1,4),c_d_11(1,4), c_d_8(1,4),c_d_5(1,4),c_d_2(1,4)];

C_d_3_true_2 = [c_d_24(1,1),c_d_21(1,1),c_d_18(1,1),c_d_15(1,1),c_d_12(1,1),...
    c_d_9(1,1),c_d_6(1,1),c_d_3(1,1),c_d_24(1,2),c_d_21(1,2),c_d_18(1,2),...
    c_d_15(1,2),c_d_12(1,2), c_d_9(1,2),c_d_6(1,2),c_d_3(1,2),...
    c_d_24(1,3),c_d_21(1,3),c_d_18(1,3),c_d_15(1,3),c_d_12(1,3),...
    c_d_9(1,3),c_d_6(1,3),c_d_3(1,3),c_d_24(1,4),c_d_21(1,4),c_d_18(1,4),...
    c_d_15(1,4),c_d_12(1,4), c_d_9(1,4),c_d_6(1,4),c_d_3(1,4)];
%% C_d plots
% naca data
naca_cd = [.01,.02,.02,.023,.038,.042,.08,.12];
naca_aoa_cd = [-8,-4,-2,0,2,4,8,12];
figure(2)
p(1) = plot(aoa,C_d_1_true_2(1,:),':ro','DisplayName','C_d at 9 m/s');
hold on
p(2) = plot(aoa,C_d_2_true_2(1,:),':bo','DisplayName','C_d at 17 m/s');
hold on 
p(3) = plot(aoa,C_d_3_true_2(1,:),':go','DisplayName','C_d at 34 m/s');
hold on
p(4) = plot(naca_aoa_cd,naca_cd,':ko','DisplayName','NACA ClarkY-14 Data');
title('$Coefficient\:of\:Drag\:Versus\:Angle\:of\:Attack\:$',...
    'Interpreter','latex')
xlabel('$Angle\:of\:Attack(degrees)$','Interpreter','latex')
ylabel('$C_d$','Interpreter','latex')
legend([p(1),p(2),p(3),p(4)])
%% Comparing slopes to NACA data
p1 = polyfit(aoa(8:16),C_l_3_true_2(8:16),1);
p2 = polyfit(naca_aoa(1:4),naca_cl(1:4),1);
percent_diff_cl = abs((p2(1)-p1(1))/(p2(1)));

p3 = polyfit(aoa(16:24),C_d_3_true_2(16:24),1);
p4 = polyfit(naca_aoa_cd(4:7),naca_cd(4:7),1);
percent_diff_cd = abs((p4(1)-p3(1))/(p4(1)));
