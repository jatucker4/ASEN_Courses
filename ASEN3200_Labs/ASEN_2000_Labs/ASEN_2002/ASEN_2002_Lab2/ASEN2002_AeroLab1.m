%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2002 Aero Lab 1
%
% Group Members
% 1) Johnathan Tucker
% 2) Adam Elsayed
% 3) Aiden Kirby
% 4) Cameron Kratt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all;
close all;
clc;
%% Fetch data from excel file Groups 1-2(Pitot Static)

[atm_pressure,atm_temp, airspeed_press, auxillary_press, eld_x, eld_y, Voltage] = read_data_xcel('VelocityVoltage_S012_G01.csv');
j = 1;
for i = 2:length(Voltage)
    if Voltage(i,1)>Voltage(i-1,1)
        index_change(j,1) = i;
        j = j+1;
    end
end
atm_pressure_F1 = mean(atm_pressure);
atm_temp_1 = mean(atm_temp);

airspeed_press_v1_f1 = mean(airspeed_press(1:index_change(1,1)-1, 1));
airspeed_press_v2_f1 = mean(airspeed_press((index_change(1,1)+1):(index_change(2,1)-1), 1));
airspeed_press_v3_f1 = mean(airspeed_press((index_change(2,1)+1):(index_change(3,1)-1), 1));
airspeed_press_v4_f1 = mean(airspeed_press((index_change(3,1)+1):(index_change(4,1)-1), 1));
airspeed_press_v5_f1 = mean(airspeed_press(index_change(4,1)+1:end, 1));

%This is for 1 volt
exit_velocity_v1_f1 = pitot_static(airspeed_press_v1_f1,atm_temp_1,atm_pressure_F1);
%This is for 3 volts
exit_velocity_v2_f1 = pitot_static(airspeed_press_v2_f1,atm_temp_1,atm_pressure_F1);
%This is for 5 volts
exit_velocity_v3_f1 = pitot_static(airspeed_press_v3_f1,atm_temp_1,atm_pressure_F1);
%This is for 7 volts
exit_velocity_v4_f1 = pitot_static(airspeed_press_v4_f1,atm_temp_1,atm_pressure_F1);
%This is for 9 volts
exit_velocity_v5_f1 = pitot_static(airspeed_press_v5_f1,atm_temp_1,atm_pressure_F1);

%% Fetching and Calculating for groups 5-6(Pitot-Static)
[atm_pressure,atm_temp, airspeed_press, auxillary_press, eld_x, eld_y, Voltage] = read_data_xcel('VelocityVoltage_S012_G05.txt');
j = 1;
for i = 2:length(Voltage)
    if Voltage(i,1)>Voltage(i-1,1)
        index_change(j,1) = i;
        j = j+1;
    end
end
atm_pressure_F3 = mean(atm_pressure);
atm_temp_3 = mean(atm_temp);

airspeed_press_v1_f3 = mean(airspeed_press(1:index_change(1,1)-1, 1));
airspeed_press_v2_f3 = mean(airspeed_press((index_change(1,1)+1):(index_change(2,1)-1), 1));
airspeed_press_v3_f3 = mean(airspeed_press((index_change(2,1)+1):(index_change(3,1)-1), 1));
airspeed_press_v4_f3 = mean(airspeed_press((index_change(3,1)+1):(index_change(4,1)-1), 1));
airspeed_press_v5_f3 = mean(airspeed_press(index_change(4,1)+1:end, 1));

%This is for 2 volt
exit_velocity_v1_f3 = pitot_static(airspeed_press_v1_f3,atm_temp_3,atm_pressure_F3);
%This is for 4 volts
exit_velocity_v2_f3 = pitot_static(airspeed_press_v2_f3,atm_temp_3,atm_pressure_F3);
%This is for 6 volts
exit_velocity_v3_f3 = pitot_static(airspeed_press_v3_f3,atm_temp_3,atm_pressure_F3);
%This is for 8 volts
exit_velocity_v4_f3 = pitot_static(airspeed_press_v4_f3,atm_temp_3,atm_pressure_F3);
%This is for 10 volts
exit_velocity_v5_f3 = pitot_static(airspeed_press_v5_f3,atm_temp_3,atm_pressure_F3);
%% Fetch Data and Calculations for Groups 9-10(Pitot-Static)
[atm_pressure,atm_temp, airspeed_press, auxillary_press, eld_x, eld_y, Voltage] = read_data_xcel('VelocityVoltage_S012_G09.csv');
j = 1;
for i = 2:length(Voltage)
    if Voltage(i,1)>Voltage(i-1,1)
        index_change(j,1) = i;
        j = j+1;
    end
end
atm_pressure_F4 = mean(atm_pressure);
atm_temp_4 = mean(atm_temp);

airspeed_press_v1_f4 = mean(airspeed_press(1:index_change(1,1)-1, 1));
airspeed_press_v2_f4 = mean(airspeed_press((index_change(1,1)+1):(index_change(2,1)-1), 1));
airspeed_press_v3_f4 = mean(airspeed_press((index_change(2,1)+1):(index_change(3,1)-1), 1));
airspeed_press_v4_f4 = mean(airspeed_press((index_change(3,1)+1):(index_change(4,1)-1), 1));
airspeed_press_v5_f4 = mean(airspeed_press(index_change(4,1)+1:end, 1));

%This is for 1.5 volt
exit_velocity_v1_f4 = pitot_static(airspeed_press_v1_f4,atm_temp_4,atm_pressure_F4);
%This is for 3.5 volts
exit_velocity_v2_f4 = pitot_static(airspeed_press_v2_f4,atm_temp_4,atm_pressure_F4);
%This is for 5.5 volts
exit_velocity_v3_f4 = pitot_static(airspeed_press_v3_f4,atm_temp_4,atm_pressure_F4);
%This is for 7.5 volts
exit_velocity_v4_f4 = pitot_static(airspeed_press_v4_f4,atm_temp_4,atm_pressure_F4);
%This is for 9.5 volts
exit_velocity_v5_f4 = pitot_static(airspeed_press_v5_f4,atm_temp_4,atm_pressure_F4);
%% Fetch Data and do Calculations for groups 13-14(Pitot-Static)
[atm_pressure,atm_temp, airspeed_press, auxillary_press, eld_x, eld_y, Voltage] = read_data_xcel('VelocityVoltage_S012_G13.csv');
j = 1;
for i = 2:length(Voltage)
    if Voltage(i,1)>Voltage(i-1,1)
        index_change(j,1) = i;
        j = j+1;
    end
end
atm_pressure_F5 = mean(atm_pressure);
atm_temp_5 = mean(atm_temp);

airspeed_press_v1_f5 = mean(airspeed_press(1:index_change(1,1)-1, 1));
airspeed_press_v2_f5 = mean(airspeed_press((index_change(1,1)+1):(index_change(2,1)-1), 1));
airspeed_press_v3_f5 = mean(airspeed_press((index_change(2,1)+1):(index_change(3,1)-1), 1));
airspeed_press_v4_f5 = mean(airspeed_press((index_change(3,1)+1):(index_change(4,1)-1), 1));
airspeed_press_v5_f5 = mean(airspeed_press(index_change(4,1)+1:end, 1));

%This is for 0.5 volt
exit_velocity_v1_f5 = pitot_static(airspeed_press_v1_f5,atm_temp_5,atm_pressure_F5);
%This is for 2.5 volts
exit_velocity_v2_f5 = pitot_static(airspeed_press_v2_f5,atm_temp_5,atm_pressure_F5);
%This is for 4.5 volts
exit_velocity_v3_f5 = pitot_static(airspeed_press_v3_f5,atm_temp_5,atm_pressure_F5);
%This is for 6.5 volts
exit_velocity_v4_f5 = pitot_static(airspeed_press_v4_f5,atm_temp_5,atm_pressure_F5);
%This is for 8.5 volts
exit_velocity_v5_f5 = pitot_static(airspeed_press_v5_f5,atm_temp_5,atm_pressure_F5);

%% Total exit velocity matrix ordered by voltage
exit_velocity_final = [exit_velocity_v1_f5; exit_velocity_v1_f1; exit_velocity_v1_f4;...
    exit_velocity_v1_f3; exit_velocity_v2_f5; exit_velocity_v2_f1; exit_velocity_v2_f4;...
    exit_velocity_v2_f3; exit_velocity_v3_f5; exit_velocity_v3_f1; exit_velocity_v3_f4;...
    exit_velocity_v3_f3; exit_velocity_v4_f5; exit_velocity_v4_f1; exit_velocity_v4_f4;...
    exit_velocity_v4_f3; exit_velocity_v5_f5; exit_velocity_v5_f1; exit_velocity_v5_f4;...
    exit_velocity_v5_f3];
%% Fetch Data and do Calculations for 3-4 (Venturi-Tube)
[atm_pressure,atm_temp, airspeed_press, auxillary_press, eld_x, eld_y, Voltage] = read_data_xcel('VelocityVoltage_S012_G03.csv');
j = 1;
for i = 2:length(Voltage)
    if Voltage(i,1)>Voltage(i-1,1)
        index_change(j,1) = i;
        j = j+1;
    end
end
atm_pressure_F6 = mean(atm_pressure);
atm_temp_6 = mean(atm_temp);

airspeed_press_v1_f6 = mean(airspeed_press(1:index_change(1,1)-1, 1));
airspeed_press_v2_f6 = mean(airspeed_press((index_change(1,1)+1):(index_change(2,1)-1), 1));
airspeed_press_v3_f6 = mean(airspeed_press((index_change(2,1)+1):(index_change(3,1)-1), 1));
airspeed_press_v4_f6 = mean(airspeed_press((index_change(3,1)+1):(index_change(4,1)-1), 1));
airspeed_press_v5_f6 = mean(airspeed_press(index_change(4,1)+1:end, 1));

%This is for 1 volt
exit_velocity_v1_f6 = venturi_tube(airspeed_press_v1_f6,atm_temp_6,atm_pressure_F6);
%This is for 3 volts
exit_velocity_v2_f6 = venturi_tube(airspeed_press_v2_f6,atm_temp_6,atm_pressure_F6);
%This is for 5 volts
exit_velocity_v3_f6 = venturi_tube(airspeed_press_v3_f6,atm_temp_6,atm_pressure_F6);
%This is for 7 volts
exit_velocity_v4_f6 = venturi_tube(airspeed_press_v4_f6,atm_temp_6,atm_pressure_F6);
%This is for 9 volts
exit_velocity_v5_f6 = venturi_tube(airspeed_press_v5_f6,atm_temp_6,atm_pressure_F6);

%% Fetch Data and do Calculations for Groups 7-8 (Venturi - Tube)
[atm_pressure,atm_temp, airspeed_press, auxillary_press, eld_x, eld_y, Voltage] = read_data_xcel('VelocityVoltage_S012_G07.csv');
j = 1;
for i = 2:length(Voltage)
    if Voltage(i,1)>Voltage(i-1,1)
        index_change(j,1) = i;
        j = j+1;
    end
end
atm_pressure_F7 = mean(atm_pressure);
atm_temp_7 = mean(atm_temp);

airspeed_press_v1_f7 = mean(airspeed_press(1:index_change(1,1)-1, 1));
airspeed_press_v2_f7 = mean(airspeed_press((index_change(1,1)+1):(index_change(2,1)-1), 1));
airspeed_press_v3_f7 = mean(airspeed_press((index_change(2,1)+1):(index_change(3,1)-1), 1));
airspeed_press_v4_f7 = mean(airspeed_press((index_change(3,1)+1):(index_change(4,1)-1), 1));
airspeed_press_v5_f7 = mean(airspeed_press(index_change(4,1)+1:end, 1));

%This is for 1 volt
exit_velocity_v1_f7 = venturi_tube(airspeed_press_v1_f7,atm_temp_7,atm_pressure_F7);
%This is for 3 volts
exit_velocity_v2_f7 = venturi_tube(airspeed_press_v2_f7,atm_temp_7,atm_pressure_F7);
%This is for 5 volts
exit_velocity_v3_f7 = venturi_tube(airspeed_press_v3_f7,atm_temp_7,atm_pressure_F7);
%This is for 7 volts
exit_velocity_v4_f7 = venturi_tube(airspeed_press_v4_f7,atm_temp_7,atm_pressure_F7);
%This is for 9 volts
exit_velocity_v5_f7 = venturi_tube(airspeed_press_v5_f7,atm_temp_7,atm_pressure_F7);
%% Fetch Data and Calculate for Groups 11-12 (Venturi - Tube)
[atm_pressure,atm_temp, airspeed_press, auxillary_press, eld_x, eld_y, Voltage] = read_data_xcel('VelocityVoltage_S012_G11.csv');
j = 1;
for i = 2:length(Voltage)
    if Voltage(i,1)>Voltage(i-1,1)
        index_change(j,1) = i;
        j = j+1;
    end
end
atm_pressure_F8 = mean(atm_pressure);
atm_temp_8 = mean(atm_temp);

airspeed_press_v1_f8 = mean(airspeed_press(1:index_change(1,1)-1, 1));
airspeed_press_v2_f8 = mean(airspeed_press((index_change(1,1)+1):(index_change(2,1)-1), 1));
airspeed_press_v3_f8 = mean(airspeed_press((index_change(2,1)+1):(index_change(3,1)-1), 1));
airspeed_press_v4_f8 = mean(airspeed_press((index_change(3,1)+1):(index_change(4,1)-1), 1));
airspeed_press_v5_f8 = mean(airspeed_press(index_change(4,1)+1:end, 1));

%This is for 1 volt
exit_velocity_v1_f8 = venturi_tube(airspeed_press_v1_f8,atm_temp_8,atm_pressure_F8);
%This is for 3 volts
exit_velocity_v2_f8 = venturi_tube(airspeed_press_v2_f8,atm_temp_8,atm_pressure_F8);
%This is for 5 volts
exit_velocity_v3_f8 = venturi_tube(airspeed_press_v3_f8,atm_temp_8,atm_pressure_F8);
%This is for 7 volts
exit_velocity_v4_f8 = venturi_tube(airspeed_press_v4_f8,atm_temp_8,atm_pressure_F8);
%This is for 9 volts
exit_velocity_v5_f8 = venturi_tube(airspeed_press_v5_f8,atm_temp_8,atm_pressure_F8);
%% Fetch Data and Do Calculations for Groups 15-16 (Venturi - Tube)
[atm_pressure,atm_temp, airspeed_press, auxillary_press, eld_x, eld_y, Voltage] = read_data_xcel('VelocityVoltage_S012_G15.csv');
j = 1;
for i = 2:length(Voltage)
    if Voltage(i,1)>Voltage(i-1,1)
        index_change(j,1) = i;
        j = j+1;
    end
end
atm_pressure_F9 = mean(atm_pressure);
atm_temp_9 = mean(atm_temp);

airspeed_press_v1_f9 = mean(airspeed_press(1:index_change(1,1)-1, 1));
airspeed_press_v2_f9 = mean(airspeed_press((index_change(1,1)+1):(index_change(2,1)-1), 1));
airspeed_press_v3_f9 = mean(airspeed_press((index_change(2,1)+1):(index_change(3,1)-1), 1));
airspeed_press_v4_f9 = mean(airspeed_press((index_change(3,1)+1):(index_change(4,1)-1), 1));
airspeed_press_v5_f9 = mean(airspeed_press(index_change(4,1)+1:end, 1));

%This is for 1 volt
exit_velocity_v1_f9 = venturi_tube(airspeed_press_v1_f9,atm_temp_9,atm_pressure_F9);
%This is for 3 volts
exit_velocity_v2_f9 = venturi_tube(airspeed_press_v2_f9,atm_temp_9,atm_pressure_F9);
%This is for 5 volts
exit_velocity_v3_f9 = venturi_tube(airspeed_press_v3_f9,atm_temp_9,atm_pressure_F9);
%This is for 7 volts
exit_velocity_v4_f9 = venturi_tube(airspeed_press_v4_f9,atm_temp_9,atm_pressure_F9);
%This is for 9 volts
exit_velocity_v5_f9 = venturi_tube(airspeed_press_v5_f9,atm_temp_9,atm_pressure_F9);

%% Final velocity vector using the airspeed pressure transducer and venturi tube calculation
exit_velocity_final_2 = [exit_velocity_v1_f9; exit_velocity_v1_f6; exit_velocity_v1_f8;...
    exit_velocity_v1_f7; exit_velocity_v2_f9; exit_velocity_v2_f6; exit_velocity_v2_f8;...
    exit_velocity_v2_f7; exit_velocity_v3_f9; exit_velocity_v3_f6; exit_velocity_v3_f8;...
    exit_velocity_v3_f7; exit_velocity_v4_f9; exit_velocity_v4_f6; exit_velocity_v4_f8;...
    exit_velocity_v4_f7; exit_velocity_v5_f9; exit_velocity_v5_f6; exit_velocity_v5_f8;...
    exit_velocity_v5_f7];
%% Boundary Layer Port 2 - 0,1,2,3,4,5,6,7,8,9,10,Center
[y_value_2(1,1),freestream1] = boundary_layer('BoundaryLayer_S012_G01.csv');

%% Boundary Layer Port 2 - 0, 0.5...9.5, Center
[y_value_2(1,2),temp] = boundary_layer('BoundaryLayer_S012_G03.csv');

%% Now for Boundary Layer Port 5 - 0,1,2,3,4,5,6,7,8,9,10,Center
[y_value_5(1,1),temp] = boundary_layer('BoundaryLayer_S011_G05.csv');

%% Now for Boundary Layer Port 5 - 0,0.5....9.5,Center
[y_value_5(1,2),temp] = boundary_layer('BoundaryLayer_s011_G07.csv');

%% Boundary Layer Port 8 - 0,1,2,3,4,5,6,7,8,9,10,Center
[y_value_8(1,1),temp] = boundary_layer('BoundaryLayer_S012_G09.csv');

%% Boundary Layer Port 8 - 0, 0.5...9.5, Center
[y_value_8(1,2),temp] = boundary_layer('BoundaryLayer_S012_G11.csv');

%% Boundary Layer Port 11 - 0,1,2,3,4,5,6,7,8,9,10,Center
[y_value_11(1,1),temp] = boundary_layer('BoundaryLayer_S011_G13.csv');

%% Boundary Layer Port 11 - 0, 0.5...9.5, Center
[y_value_11(1,2),freestream2] = boundary_layer('BoundaryLayer_S012_G15.csv');

%% Water manometer
Pdiff = [0.2, 0.6, 2.05, 3.8, 5.8];
for i = 1:length(Pdiff)
V_manometer(1,i) = venturi_tube((Pdiff(1,i)*284.84),atm_temp_3,atm_pressure_F3);
end

%% Calculate V2 with boundary layer influence using continuity eqn
V1 = freestream1;
A_1 = (12*12)/1550.003;
A_2 = 0;

%% Plots for pitot static
voltages = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10];
voltages_2 = [2,4,6,8,10];
p = polyfit(voltages, exit_velocity_final', 1);
p2 = polyfit(voltages, exit_velocity_final_2', 1);
p3 = polyfit(voltages_2, V_manometer, 1);
y = p(1,1)*voltages + p(1,2);
y_2 = p2(1,1)*voltages + p2(1,2);
plot(voltages, exit_velocity_final)
hold on
% plot(voltages_2, V_manometer)
plot(voltages, y)
% hold on
% plot(voltages, exit_velocity_final_2)
xlabel('Voltages (V)')
ylabel('Airspeed (m/s)')
title('Airspeed versus Voltage Graph')
legend('Pitot-Static Calculation', 'Water Manometer', 'Venturi-Tube Calculation')
