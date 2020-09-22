function [y_value,v_inf]  = boundary_layer(inputfile)
[atm_pressure,atm_temp, airspeed_press, auxillary_press, eld_x, eld_y, Voltage] = read_data_xcel(inputfile);
j = 1;
index_change = [501, 1001, 1501, 2001, 2501, 3001, 3501, 4001, 4501, 5001, 5501];

atm_pressure_F1 = mean(atm_pressure);
atm_temp_1 = mean(atm_temp);

airspeed_press_1(1,1) = mean(airspeed_press(1:index_change(1,1)-1, 1));
airspeed_press_1(1,2) = mean(airspeed_press((index_change(1,1)+1):(index_change(1,2)), 1));
airspeed_press_1(1,3) = mean(airspeed_press((index_change(1,2)+1):(index_change(1,3)), 1));
airspeed_press_1(1,4) = mean(airspeed_press((index_change(1,3)+1):(index_change(1,4)), 1));
airspeed_press_1(1,5) = mean(airspeed_press(index_change(1,4)+1:(index_change(1,5)), 1));
airspeed_press_1(1,6) = mean(airspeed_press(index_change(1,5)+1:(index_change(1,6)), 1));
airspeed_press_1(1,7) = mean(airspeed_press(index_change(1,6)+1:(index_change(1,7)), 1));
airspeed_press_1(1,8) = mean(airspeed_press(index_change(1,7)+1:(index_change(1,8)), 1));
airspeed_press_1(1,9) = mean(airspeed_press(index_change(1,8)+1:(index_change(1,9)), 1));
airspeed_press_1(1,10) = mean(airspeed_press(index_change(1,9)+1:(index_change(1,10)), 1));
airspeed_press_1(1,11) = mean(airspeed_press(index_change(1,10)+1:(index_change(1,11)), 1));
airspeed_press_1(1,12) = mean(airspeed_press(index_change(1,11)+1:end, 1));

% Calculate free stream velocity
v_inf = pitot_static(airspeed_press_1(1,12), atm_temp_1, atm_pressure_F1);

for i = 1:length(airspeed_press_1)-1
    V(1,i) = pitot_static(airspeed_press_1(1,i), atm_temp_1, atm_pressure_F1);
end

V = V - v_inf*.95;
boundary_index = find(V == min(V),1);
if boundary_index == 1
    y_value = mean(eld_y(1:501,1));
else 
    y_value = mean(eld_y(index_change(boundary_index-1) : index_change(boundary_index)));
end
end

